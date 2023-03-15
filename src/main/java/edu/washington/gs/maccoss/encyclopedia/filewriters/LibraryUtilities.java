package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryEntryCleaner;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Quadruplet;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.LibraryEntryModifier;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakChromatogram;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
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
	
	public static void extractSampleSpecificLibraries(ProgressIndicator progress, File saveDir, Optional<HashMap<String, double[]>> targetTransitions, LibraryInterface library, final SearchParameters parameters) throws IOException, SQLException, DataFormatException {
		List<Path> sources=library.getSourceFiles();
		saveDir.mkdirs();
		for (Path path : sources) {
			String name=path.getFileName().toString();
			progress.update("Extracting out spectra associated with ["+name+"]...");
			Logger.logLine("Extracting out spectra associated with ["+name+"]...");
			File saveFile=new File(saveDir, name+LibraryFile.DLIB);
			extractSampleSpecificLibrary(progress, saveFile, name, true, targetTransitions.isPresent()?targetTransitions.get():null, library, parameters);
		}
		progress.update("Finished extracting "+sources.size()+" batch-specific libraries!");
	}
		
	private static void extractSampleSpecificLibrary(final ProgressIndicator progress, final File saveFile, String sourceFile, boolean useBestFragmentation, HashMap<String, double[]> targetTransitions, LibraryInterface library, final SearchParameters parameters) throws IOException, SQLException, DataFormatException {
		LibraryFile saveLibrary=new LibraryFile();
		saveLibrary.openFile();
		MassTolerance tolerance=parameters.getFragmentTolerance();
		
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
			LibraryEntry thisSamplesEntry=entry.getValue();
			LibraryEntry bestEntry=bestSource.get(key);
			String peptideModSeq=thisSamplesEntry.getPeptideModSeq();
			
			double[] targets=null; // targets stays null if no transitions are targetted
			if (targetTransitions!=null) {
				targets=targetTransitions.get(peptideModSeq);
				if (targets==null) {
					Logger.errorLine("Missing transitions for "+peptideModSeq+", skipping peptide");
					continue; // ignore peptides that don't have transitions
				} else if (targets.length==0 ) {
					Logger.errorLine("Zero targeted transitions for "+peptideModSeq+", skipping peptide");
					continue; // ignore peptides that don't have transitions
				}
			}

			double[] massArray;
			float[] intensityArray;
			float[] correlationArray;
			boolean[] quantifiedIonsArray;
			if (bestEntry!=null&&bestEntry.getScore()<thisSamplesEntry.getScore()) {
				massArray = bestEntry.getMassArray();
				intensityArray = bestEntry.getIntensityArray();
				correlationArray = bestEntry.getCorrelationArray();
				quantifiedIonsArray = bestEntry.getQuantifiedIonsArray();
			} else {
				massArray = thisSamplesEntry.getMassArray();
				intensityArray = thisSamplesEntry.getIntensityArray();
				correlationArray = thisSamplesEntry.getCorrelationArray();
				quantifiedIonsArray = thisSamplesEntry.getQuantifiedIonsArray();
			}
			
			if (targets!=null) {
				// reset quantifiedIonsArray for targets
				quantifiedIonsArray=new boolean[massArray.length];
				TDoubleArrayList missing=new TDoubleArrayList();
				for (int i = 0; i < targets.length; i++) {
					Optional<Integer> index=tolerance.getIndex(massArray, targets[i]);
					if (index.isPresent()) {
						quantifiedIonsArray[index.get()]=true;
					} else {
						missing.add(targets[i]);
					}
				}
				
				if (missing.size()>0) {
					// if targets are missing, then make sure they're added
					ArrayList<PeakChromatogram> peaks=new ArrayList<PeakChromatogram>();
					for (int i = 0; i < massArray.length; i++) {
						peaks.add(new PeakChromatogram(massArray[i], intensityArray[i], correlationArray[i], quantifiedIonsArray[i]));
					}
					
					for (double mass : missing.toArray()) {
						peaks.add(new PeakChromatogram(mass, Float.MIN_NORMAL, 0.0f, true));
					}
					
					Collections.sort(peaks);
					Quadruplet<double[], float[], float[], boolean[]> arrays=PeakChromatogram.toChromatogramArrays(peaks);
					
					massArray=arrays.x;
					intensityArray=arrays.y;
					correlationArray=arrays.z;
					quantifiedIonsArray=arrays.w;
				}
			}

			toWrite.add(thisSamplesEntry.updateMS2(massArray, intensityArray, correlationArray, quantifiedIonsArray));
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
	
	public static void subsetLibrary(final File saveFile, final float rtMinSec, final float rtMaxSec, final double mzMin, final double mzMax, final HashSet<String> targets, LibraryInterface library)
			throws IOException, SQLException, DataFormatException {
		LibraryFile saveLibrary=new LibraryFile();
		saveLibrary.openFile();
		
		ArrayList<LibraryEntry> toWrite=new ArrayList<>();
		for (LibraryEntry entry : library.getAllEntries(false, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()))) {
			if (rtMinSec<=entry.getRetentionTime()&&rtMaxSec>=entry.getRetentionTime()&&mzMin<=entry.getPrecursorMZ()&&mzMax>=entry.getPrecursorMZ()) {
				if (matchesKeyword(targets, entry)) {
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
	
	private static boolean matchesKeyword(HashSet<String> targets, LibraryEntry entry) {
		if (targets.size()==0) return true;
		if (targets.contains(entry.getPeptideSeq())) return true;
		if (targets.contains(entry.getPeptideModSeq())) return true;
		if (targets.contains(entry.getSource())) return true;
		
		for (String accession : entry.getAccessions()) {
			if (targets.contains(accession)) return true;
		}
		return false;
	}
	
	public static LibraryFile mergeLibraries(ProgressIndicator progress, ArrayList<File> files, File saveFile, boolean rtAlign, boolean removeDuplicates, boolean higherScoresAreBetter) throws IOException, SQLException, DataFormatException {
		HashMap<String, ArrayList<LibraryEntry>> groupedEntries=new HashMap<>();
		int totalEntries=0;
		int count=0;
		for (File elibFile : files) {
			count++;
			LibraryFile library=new LibraryFile();
			library.openFile(elibFile);
			ArrayList<LibraryEntry> localEntries = library.getAllEntries(false,  new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()));
			//groupedEntries.put(elibFile.getName(), localEntries);
			for (LibraryEntry entry : localEntries) {
				ArrayList<LibraryEntry> list=groupedEntries.get(entry.getSource());
				if (list==null) {
					list=new ArrayList<>();
					groupedEntries.put(entry.getSource(), list);
				}
				list.add(entry);
				totalEntries++;
			}
			progress.update("Found "+localEntries.size()+" entries from "+elibFile.getName(), count/(files.size()+1.0f));
			Logger.logLine("Found "+localEntries.size()+" entries from "+elibFile.getName()+", "+totalEntries+" total entries from "+groupedEntries.size()+" sources...");
			library.close();
		}
		
		ArrayList<LibraryEntry> allEntries;
		if (rtAlign) {
			progress.update("Correcting retention times...");
			allEntries=LibraryEntryCleaner.correctRTs(groupedEntries, saveFile);
		} else {
			allEntries=new ArrayList<>();
			for (ArrayList<LibraryEntry> list : groupedEntries.values()) {
				allEntries.addAll(list);
			}
		}
		if (removeDuplicates) {
			allEntries=LibraryEntryCleaner.removeDuplicateEntries(allEntries, higherScoresAreBetter);
			progress.update("Removing duplicates...");
		}

		LibraryFile saveLibrary=new LibraryFile();
		saveLibrary.openFile();
		saveLibrary.dropIndices();
		saveLibrary.addEntries(allEntries);
		saveLibrary.addProteinsFromEntries(allEntries);
		saveLibrary.createIndices();
		saveLibrary.saveAsFile(saveFile);
		
		saveLibrary.close();
		progress.update("Saved "+saveFile.getName()+", "+allEntries.size()+" total", 1.0f);
		Logger.logLine("Saved "+saveFile.getName()+", "+allEntries.size()+" total");
		return saveLibrary;
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
