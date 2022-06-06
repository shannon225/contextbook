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
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

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

	/**
	 * @return A <emph>closed</emph> {@code LibraryFile} instance pointing to {@code elibFile}
	 *         containing the results of conversion.
	 */
	public static LibraryFile convertFromSpectronautCSV(File csvFile, File fastaFile, SearchParameters parameters) {
		String absolutePath=csvFile.getAbsolutePath();
		File libraryFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.DLIB);
		return convertFromSpectronautCSV(csvFile, fastaFile, libraryFile, parameters);
	}

	/**
	 * @return A <emph>closed</emph> {@code LibraryFile} instance pointing to {@code elibFile}
	 *         containing the results of conversion.
	 */
	public static LibraryFile convertFromSpectronautCSV(File csvFile, File fastaFile, File libraryFile, SearchParameters parameters) {
		AminoAcidConstants aaConstants=parameters.getAAConstants();
		try {
			String sourceFile = csvFile.getName();
			Logger.logLine("Started parsing entries in "+sourceFile+"...");
			
			final ArrayList<ImmutablePeptideEntry> peptides=new ArrayList<ImmutablePeptideEntry>();
			TableParserMuscle muscle=new TableParserMuscle() {
				private PeptideEntry lastPeptide=null;
				private String lastGroup=null;
				@Override
				public void processRow(Map<String, String> row) {
					try {
						String peptideModSeq=OpenSwathTSVToLibraryConverter.getFromMap(row, "ModifiedPeptide");
						String chargeString=OpenSwathTSVToLibraryConverter.getFromMap(row, "PrecursorCharge");
						double productMz=Double.parseDouble(OpenSwathTSVToLibraryConverter.getFromMap(row, "FragmentMz"));
						float libraryIntensity=Float.parseFloat(OpenSwathTSVToLibraryConverter.getFromMap(row, "RelativeIntensity", "RelativeFragmentIntensity"));
	
						String group=peptideModSeq+"_"+chargeString+"H";
						
						if (!group.equals(lastGroup)) {
							byte charge=Byte.parseByte(chargeString);
							float iRT=Float.parseFloat(OpenSwathTSVToLibraryConverter.getFromMap(row, "iRT"));
							
							if (lastPeptide!=null) peptides.add(new ImmutablePeptideEntry(lastPeptide));
							
							lastPeptide=new PeptideEntry(parseMods(peptideModSeq), charge, iRT*60f, sourceFile);
							lastGroup=group;
							
							if (peptides.size()%10000==0) {
								Logger.logLine("Read "+peptides.size()+" entries...");
							}
						}
						lastPeptide.addPeak(new Peak(productMz, libraryIntensity));

					} catch (Exception e) {
						Logger.errorLine("Error parsing Spectronaut CSV:");
						Logger.errorException(e);
						Logger.errorLine("Spectronaut CSV parsing requires the following columns:\n" +
								" 1) ModifiedPeptide: ["+OpenSwathTSVToLibraryConverter.getFromMap(row, "ModifiedPeptide")+"]\n" + 
								" 2) PrecursorCharge: ["+OpenSwathTSVToLibraryConverter.getFromMap(row, "PrecursorCharge")+"]\n" + 
								" 3) FragmentMz: ["+OpenSwathTSVToLibraryConverter.getFromMap(row, "FragmentMz")+"]\n" + 
								" 4) RelativeIntensity: ["+OpenSwathTSVToLibraryConverter.getFromMap(row, "RelativeIntensity")+"]\n" + 
								" 5) iRT: ["+OpenSwathTSVToLibraryConverter.getFromMap(row, "iRT")+"]");
						Logger.errorLine("Other headers: "+General.toString(new ArrayList<String>(row.keySet())));
						throw new EncyclopediaException(e);
					}
				}
				
				@Override
				public void cleanup() {
					if (lastPeptide!=null) peptides.add(new ImmutablePeptideEntry(lastPeptide));
				}
			};

			String lower = sourceFile.toLowerCase();
			boolean isTSV=lower.endsWith(".xls")||lower.endsWith(".txt")||lower.endsWith(".tsv");
			if (isTSV) {
				TableParser.parseTSV(csvFile, muscle);
			} else {
				// otherwise (e.g. *.spectronaut, etc, parse as CSV)
				TableParser.parseCSV(csvFile, muscle);
			}

			Logger.logLine("Finished parsing, now processing entries...");
			return OpenSwathTSVToLibraryConverter.processPeptideEntries(sourceFile, fastaFile, libraryFile, parameters, aaConstants, peptides);

		} catch (Exception e) {
			Logger.errorLine("Error parsing Spectronaut CSV:");
			Logger.errorException(e);
			throw new EncyclopediaException(e);
		}
	}
}
