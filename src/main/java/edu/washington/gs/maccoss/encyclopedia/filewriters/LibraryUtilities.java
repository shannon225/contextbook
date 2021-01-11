package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.zip.DataFormatException;

import org.jfree.base.Library;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.LibraryEntryModifier;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class LibraryUtilities {
	public static void modifyLibrary(final File saveFile, TCharDoubleHashMap modMasses, boolean isFixed, LibraryInterface library) throws IOException, SQLException, DataFormatException {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		
		LibraryFile saveLibrary=new LibraryFile();
		saveLibrary.openFile();
		
		ArrayList<LibraryEntry> toWrite=new ArrayList<>();
		int count=0;
		for (LibraryEntry entry : library.getAllEntries(false, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()))) {
			count++;
			if (isFixed) {
				toWrite.add(LibraryEntryModifier.modifyModelAtEverySite(entry, modMasses, false, parameters));
			} else {
				toWrite.addAll(LibraryEntryModifier.modifyModelAtEachSite(entry, modMasses, true, parameters));
			}
		}
		Logger.logLine("Created "+toWrite.size()+" peptides from "+count+" target sequences. Writing to ["+saveFile.getAbsolutePath()+"]...");
		
		saveLibrary.dropIndices();
		saveLibrary.addEntries(toWrite);
		saveLibrary.addProteinsFromEntries(toWrite);
		saveLibrary.createIndices();
		saveLibrary.saveAsFile(saveFile);
		
		saveLibrary.close();
	}
	
	public static void extractSampleSpecificLibrary(final File saveFile, String sourceFile, boolean useBestFragmentation, LibraryInterface library) throws IOException, SQLException, DataFormatException {
		LibraryFile saveLibrary=new LibraryFile();
		saveLibrary.openFile();
		
		HashMap<String, LibraryEntry> targetSource=new HashMap<>();
		HashMap<String, LibraryEntry> bestSource=new HashMap<>();		
		for (LibraryEntry entry : library.getAllEntries(false, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()))) {
			if (entry.getRetentionTime()>=0) {
				String key=entry.getPeptideModSeq()+"+"+entry.getPrecursorCharge();
				if (sourceFile.equals(entry.getSource())) {
					updateEntryIfBetter(targetSource, entry, key);
				}
				if (useBestFragmentation) {
					updateEntryIfBetter(bestSource, entry, key);
				}
			}
		}
		
		ArrayList<LibraryEntry> toWrite=new ArrayList<>();
		for (Entry<String, LibraryEntry> entry : targetSource.entrySet()) {
			String key=entry.getKey();
			LibraryEntry targetEntry=entry.getValue();
			LibraryEntry bestEntry=bestSource.get(key);
			
			if (bestEntry!=null&&bestEntry.getScore()<targetEntry.getScore()) {
				toWrite.add(targetEntry.updateMS2(bestEntry.getMassArray(), bestEntry.getIntensityArray(), bestEntry.getCorrelationArray()));
			} else {
				toWrite.add(targetEntry);
			}
		}
		Logger.logLine("Found "+toWrite.size()+" peptides from "+sourceFile+". Writing to ["+saveFile.getAbsolutePath()+"]...");
		
		saveLibrary.dropIndices();
		saveLibrary.addEntries(toWrite);
		saveLibrary.addProteinsFromEntries(toWrite);
		saveLibrary.createIndices();
		saveLibrary.saveAsFile(saveFile);
		
		saveLibrary.close();
	}

	private static void updateEntryIfBetter(HashMap<String, LibraryEntry> targetSource, LibraryEntry entry, String key) {
		LibraryEntry alt=targetSource.get(key);
		if (alt==null||alt.getScore()>entry.getScore()) {
			targetSource.put(key, entry);
		}
	}
	
	public static void subsetLibrary(final File saveFile, final float rtMinSec, final float rtMaxSec, final HashSet<String> targets, LibraryInterface library)
			throws IOException, SQLException, DataFormatException {
		LibraryFile saveLibrary=new LibraryFile();
		saveLibrary.openFile();
		
		ArrayList<LibraryEntry> toWrite=new ArrayList<>();
		for (LibraryEntry entry : library.getAllEntries(false, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()))) {
			if (rtMinSec<=entry.getRetentionTime()&&rtMaxSec>=entry.getRetentionTime()) {
				if (targets.size()==0||(targets.contains(entry.getPeptideSeq())||targets.contains(entry.getPeptideModSeq()))) {
					toWrite.add(entry);
				}
			}
		}
		Logger.logLine("Found "+toWrite.size()+" peptides from "+targets.size()+" target sequences. Writing to ["+saveFile.getAbsolutePath()+"]...");
		
		saveLibrary.dropIndices();
		saveLibrary.addEntries(toWrite);
		saveLibrary.addProteinsFromEntries(toWrite);
		saveLibrary.createIndices();
		saveLibrary.saveAsFile(saveFile);
		
		saveLibrary.close();
	}
	
	/*
	 * FOR MSPLIT-DIA
	 */
	public static void libraryToMGF(final File saveFile, LibraryInterface library)
			throws IOException, SQLException, DataFormatException {
		PrintWriter writer=new PrintWriter(saveFile, "UTF-8");
		
		int count=0;
		for (LibraryEntry entry : library.getAllEntries(false, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()))) {
			count++;
			writer.println("BEGIN IONS");
			writer.println("Title=Scan Number: "+entry.getSpectrumIndex()+" Retention Time: PT"+entry.getScanStartTime()+"S PROTEIN:"+PSMData.accessionsToString(entry.getAccessions()));
			writer.println("CHARGE="+entry.getPrecursorCharge());
			writer.println("PEPMASS="+(entry.getPrecursorMZ()*entry.getPrecursorCharge()-entry.getPrecursorCharge()*MassConstants.protonMass));
			writer.println("SEQ="+entry.getPeptideModSeq().replace("[", "").replace("]", ""));
			writer.println("RTINSECONDS="+entry.getScanStartTime());
			writer.println("SCANS="+Math.round(entry.getScanStartTime())*1000); // assumes scans happen less often than every 1ms

			double[] masses=entry.getMassArray();
			float[] intensities=entry.getIntensityArray();
			for (int i = 0; i < intensities.length; i++) {
				writer.println(masses[i]+" "+intensities[i]);
			}
			
			writer.println("END IONS");
			writer.println();
		}
		writer.close();
		Logger.logLine("Found "+count+" peptides. Writing to ["+saveFile.getAbsolutePath()+"]...");
	}
	
	public static void main(String[] args) throws Exception {
		File libFile=new File("/Volumes/searle_ssd/malaria/novo_yeast/libraries/uniprot_yeast_25jan2019.fasta.z2_nce33.dlib");
		File mgfFile=new File("/Volumes/searle_ssd/malaria/novo_yeast/libraries/uniprot_yeast_25jan2019.fasta.z2_nce33.mgf");
		
		LibraryFile library=new LibraryFile();
		library.openFile(libFile);
		
		libraryToMGF(mgfFile, library);
	}
}
