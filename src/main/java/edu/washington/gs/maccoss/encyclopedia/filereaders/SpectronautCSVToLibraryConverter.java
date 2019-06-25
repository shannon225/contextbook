package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.ArrayList;
import java.util.Map;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PTMMap.PostTranslationalModification;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;

public class SpectronautCSVToLibraryConverter {

	// example: _YNDLEKHVC[Carbamidomethyl (C)]SIK_
	static String parseMods(String structuredSequence) {
		if (structuredSequence.charAt(0)=='_'&&structuredSequence.charAt(structuredSequence.length()-1)=='_') {
			structuredSequence=structuredSequence.substring(1, structuredSequence.length()-1);
		}
		
		char[] ca=structuredSequence.toCharArray();
		
		ArrayList<String> aas=new ArrayList<String>();
		for (int i = 0; i < ca.length; i++) {
			if (ca[i]=='[') {
				StringBuilder sb=new StringBuilder();
				i++;
				while (ca[i]!=']') {
					sb.append(ca[i]);
					i++;
				}
				if (aas.size()==0) {
					// handling of n-termini mods assumes you can't have multiple []s in a row
					i++;
					aas.add(Character.toString(ca[i]));
				}
				String ptmDescription=sb.toString();
				String ptmName=ptmDescription;
				String ptmSpecificity;
				if (ptmDescription.indexOf(' ')>0) {
					ptmName=ptmDescription.substring(0, ptmDescription.indexOf(' '));
					ptmSpecificity=ptmDescription.substring(ptmDescription.indexOf(' ')+1);
					if (ptmSpecificity.charAt(0)=='('&&ptmSpecificity.charAt(ptmSpecificity.length()-1)==')') {
						ptmSpecificity=ptmSpecificity.substring(1, ptmSpecificity.length()-1);
					}
				} else {
					ptmName=ptmDescription;
					ptmSpecificity="";
				}
				
				PostTranslationalModification ptm=PTMMap.getPTM(ptmName, ptmSpecificity);
				double modificationMass = ptm.getDeltaMass();
				
				String aaString=aas.get(aas.size()-1);
				aas.set(aas.size()-1, aaString+(modificationMass>=0?"[+":"[")+modificationMass+"]");
			} else {
				aas.add(Character.toString(ca[i]));
			}
		}
		
		StringBuilder sb=new StringBuilder();
		for (String aa : aas) {
			sb.append(aa);
		}
		return sb.toString();
	}
	
	public static LibraryInterface convertFromSpectronautCSV(File csvFile, File fastaFile, SearchParameters parameters) {
		String absolutePath=csvFile.getAbsolutePath();
		File libraryFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.DLIB);
		return convertFromSpectronautCSV(csvFile, fastaFile, libraryFile, parameters);
	}

	//FragmentCharge,FragmentNumber,FragmentType,PrecursorCharge,PrecursorMz,RelativeIntensity,iRT,ModifiedPeptide,StrippedPeptide,FragmentLossType,FragmentMz,LabeledPeptide
	public static LibraryFile convertFromSpectronautCSV(File csvFile, File fastaFile, File libraryFile, SearchParameters parameters) {
		AminoAcidConstants aaConstants=parameters.getAAConstants();
		try {
			Logger.logLine("Started parsing entries in "+csvFile.getName()+"...");
			
			final ArrayList<ImmutablePeptideEntry> peptides=new ArrayList<ImmutablePeptideEntry>();
			TableParserMuscle muscle=new TableParserMuscle() {
				private PeptideEntry lastPeptide=null;
				private String lastGroup=null;
				@Override
				public void processRow(Map<String, String> row) {
					String peptideModSeq=OpenSwathTSVToLibraryConverter.getFromMap(row, "ModifiedPeptide");
					String chargeString=OpenSwathTSVToLibraryConverter.getFromMap(row, "PrecursorCharge");
					double productMz=Double.parseDouble(OpenSwathTSVToLibraryConverter.getFromMap(row, "FragmentMz"));
					float libraryIntensity=Float.parseFloat(OpenSwathTSVToLibraryConverter.getFromMap(row, "RelativeIntensity"));

					String group=peptideModSeq+"_"+chargeString+"H";
					
					if (!group.equals(lastGroup)) {
						byte charge=Byte.parseByte(chargeString);
						float iRT=Float.parseFloat(OpenSwathTSVToLibraryConverter.getFromMap(row, "iRT"));
						
						if (lastPeptide!=null) peptides.add(new ImmutablePeptideEntry(lastPeptide));
						
						lastPeptide=new PeptideEntry(parseMods(peptideModSeq), charge, iRT);
						lastGroup=group;
						
						if (peptides.size()%10000==0) {
							Logger.logLine("Read "+peptides.size()+" entries...");
						}
					}
					lastPeptide.addPeak(new Peak(productMz, libraryIntensity));

				}
				
				@Override
				public void cleanup() {
					if (lastPeptide!=null) peptides.add(new ImmutablePeptideEntry(lastPeptide));
				}
			};
			
			TableParser.parseCSV(csvFile, muscle);

			Logger.logLine("Finished parsing, now processing entries...");
			return OpenSwathTSVToLibraryConverter.processPeptideEntries(csvFile.getName(), fastaFile, libraryFile, parameters, aaConstants, peptides);

		} catch (Exception e) {
			Logger.errorLine("Error parsing Spectronaut CSV:");
			Logger.errorException(e);
			throw new EncyclopediaException(e);
		}
	}
}
