package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeptideXYPoint;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.TwoDimensionalKDE;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideAccessionMatchingTrie;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Quadruplet;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakChromatogram;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.IndexedObject;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import gnu.trove.list.array.TFloatArrayList;
import org.apache.commons.io.FilenameUtils;

import java.io.File;
import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.zip.DataFormatException;

public class LibraryEntryCleaner {
	public static void main(String[] args) throws Exception {
//		File originalLibraryFile=new File("/Users/searleb/Documents/iarpa/new_llnl_individual/IARPA_3clib_plus_llnl.trypsin.z3_nce33.dlib");
//		File newLibraryFile=new File("/Users/searleb/Documents/iarpa/new_llnl_individual/uniq_IARPA_3clib_plus_llnl.trypsin.z3_nce33.dlib");
//		File fastaFile=new File("/Users/searleb/Documents/iarpa/new_llnl_individual/IARPA_20201215_plus_ref_vars.fasta");
//		
//		cleanLibrary(true, originalLibraryFile, newLibraryFile, fastaFile, SearchParameterParser.getDefaultParametersObject());
		

		//File fileName=new File("/Users/searleb/Documents/grants/phosphotau mutagenesis/barcodes_prosit.dlib");
		//File newFileName=new File("/Users/searleb/Documents/grants/phosphotau mutagenesis/ecoli_barcodes_prosit.dlib");
		//File fasta=new File("/Users/searleb/Documents/grants/phosphotau mutagenesis/uniprot_ecoli_27dec2020.fasta");

		//fileName = new File("/Users/searleb/Documents/cobbs/2020dec30_cobbs/2020dec03_cobbs_cmv_inf_gpfdia_clib.elib");
		//newFileName = new File("/Users/searleb/Documents/cobbs/2020dec30_cobbs/2020dec03_cobbs_cmv_inf_gpfdia_clib_hcmv_only.elib");
		//fasta=new File("/Users/searleb/Documents/cobbs/2020dec30_cobbs/merlin.fasta");
		
		//File fileName=new File("/Users/searleb/Documents/iarpa/new_llnl_individual_2/unclean_IARPA_20210401_3clib_plus_llnl.trypsin.z3_nce33.dlib");
		//File newFileName=new File("/Users/searleb/Documents/iarpa/new_llnl_individual_2/clean_IARPA_20210401_3clib_plus_llnl.trypsin.z3_nce33.dlib");
		//File fasta=new File("/Users/searleb/Documents/iarpa/new_llnl_individual_2/IARPA_20210401.fasta_plus_llnl_ref_vars.fasta");
		
//		File fileName=new File("/Users/searleb/Documents/iarpa/IARPA_bone/10p/unclean_bone-refs-vars-Pplus16-20210615.trypsin.z3_nce33.prosit.dlib");
//		File newFileName=new File("/Users/searleb/Documents/iarpa/IARPA_bone/10p/clean_bone-refs-vars-Pplus16-20210615.trypsin.z3_nce33.prosit.dlib");
//		File fasta=new File("/Users/searleb/Documents/iarpa/IARPA_bone/10p/bone-refs-vars-20210615.combined.fasta");

		File fileName=new File("/Users/searleb/Documents/phospho/phosphopedia_hcd_combined.dlib");
		File newFileName=new File("/Users/searleb/Documents/phospho/phosphopedia_hcd_combined_cleaned.dlib");
		
		fileName=new File("/Users/searleb/Documents/OSU/damien/hela_raws/HeLa.elib");
		newFileName=new File("/Users/searleb/Documents/OSU/damien/hela_raws/Villen_HeLa_cleaned.elib");
		File fasta=new File("/Users/searleb/Documents/phospho/uniprot_human_25apr2019.fasta");

//		File fileName=new File("/Users/searleb/Documents/iarpa/sigsci/unclean_IARPA_20210401_3clib_plus_llnl_and_sigsci.trypsin.z3_nce33.dlib");
//		File newFileName=new File("/Users/searleb/Documents/iarpa/sigsci/clean_IARPA_20210401_3clib_plus_llnl_and_sigsci.trypsin.z3_nce33.dlib");
//		File fasta=new File("/Users/searleb/Documents/iarpa/sigsci/IARPA_20210401.fasta_plus_llnl_and_sigsci_ref_vars.fasta");
		
		cleanLibrary(true, true, fileName, newFileName, fasta, SearchParameterParser.getDefaultParametersObject());
	}
	
	public static boolean doesLibraryNeedCleaning(LibraryFile libraryFile) throws IOException, SQLException {
		Connection c = libraryFile.getConnection();
		try {
			Statement s=c.createStatement();
			try {
				Logger.logLine("Getting count of protein accessions...");
				ResultSet rs=s.executeQuery("select count(ProteinAccession) from peptidetoprotein");
				int count=rs.getInt(0);
				Logger.logLine("...found "+count);
				
				if (count>1) {
					Logger.logLine(libraryFile.getName()+" does not need to be cleaned!");
					return false;
				}
				rs.close();

			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
		Logger.logLine(libraryFile.getName()+" needs to be cleaned!");
		return true;
	}
	
	public static LibraryFile cleanLibrary(boolean smallerScoresAreGood, boolean respectEnzyme, File originalLibraryFile, File newLibraryFile, File fastaFile, SearchParameters parameters) throws IOException, DataFormatException, SQLException {
		LibraryFile oldLibrary=new LibraryFile();
		oldLibrary.openFile(originalLibraryFile);
		ArrayList<LibraryEntry> originalEntries=oldLibrary.getAllEntries(false, parameters.getAAConstants());
		
		HashMap<String, LibraryEntry> uniqueEntryMap=new HashMap<>();
		for (LibraryEntry entry : originalEntries) {
			entry.getAccessions().clear(); // clear accessions so they can be re-populated
			String key=entry.getPeptideModSeq()+"_"+entry.getPrecursorCharge();
			
			LibraryEntry prev=uniqueEntryMap.get(key);
			if (prev==null||(smallerScoresAreGood?(prev.getScore()>entry.getScore()):(prev.getScore()<entry.getScore()))) {
				uniqueEntryMap.put(key, entry);
			}
		}
		ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>(uniqueEntryMap.values());
		
		Logger.logLine("Found "+entries.size()+" unique peptide entries");

		if (fastaFile!=null) {
			Logger.logLine("Reading Fasta file "+fastaFile.getName());
			ArrayList<FastaEntryInterface> proteins=FastaReader.readFasta(fastaFile, parameters);
		
			Logger.logLine("Constructing trie from library peptides");
			DigestionEnzyme enzyme = respectEnzyme?parameters.getEnzyme():null;
			PeptideAccessionMatchingTrie trie=new PeptideAccessionMatchingTrie(entries, Optional.ofNullable(enzyme));
			trie.addFasta(proteins);
		}

		int[] counts=new int[21];
		for (LibraryEntry entry : entries) {
			int size=Math.min(counts.length-1, entry.getAccessions().size());
			counts[size]++;
			if (size==0) {
				System.out.println(entry.getPeptideSeq());
			}
		}
		Logger.logLine("Accession count histogram: ");
		for (int i=0; i<counts.length; i++) {
			Logger.logLine(i+" Acc\t"+counts[i]+" Counts");
		}

		if (counts[0]>0) {
			Logger.errorLine(counts[0]+" library entries can't be linked to proteins! These entries will be dropped.");
		}
		
		LibraryFile newLibrary=new LibraryFile();
		newLibrary.openFile();
		Logger.logLine("Writing library file "+newLibrary.getName());
		newLibrary.dropIndices();
		newLibrary.addEntries(entries);
		newLibrary.addProteinsFromEntries(entries);
		newLibrary.createIndices();
		newLibrary.saveAsFile(newLibraryFile);
		return newLibrary;
	}

	public static ArrayList<LibraryEntry> removeDuplicateEntries(ArrayList<LibraryEntry> originalEntries, boolean higherScoresAreBetter) {
		HashMap<String, LibraryEntry> entriesByKey=new HashMap<>();
		for (LibraryEntry entry : originalEntries) {
			String key=entry.getPeptideModSeq()+"+"+entry.getPrecursorCharge();
			LibraryEntry alt=entriesByKey.get(key);
			if (alt==null) {
				entriesByKey.put(key, entry);
			} else if (higherScoresAreBetter&&entry.getScore()>alt.getScore()) {
				entriesByKey.put(key, entry);
			} else if ((!higherScoresAreBetter)&&entry.getScore()<alt.getScore()) {
				entriesByKey.put(key, entry);
			}
		}
		ArrayList<LibraryEntry> retEntries=new ArrayList<>(entriesByKey.values());
		Collections.sort(retEntries);
		return retEntries;
	}
	
	public static ArrayList<LibraryEntry> correctRTsVsSingleSource(HashMap<String, ArrayList<LibraryEntry>> originalEntriesBySourceFile, ArrayList<LibraryEntry> singleSource, File libraryFile) {
		Logger.logLine("Comparing alignment with "+singleSource.size()+" data points");

		HashMap<String, ArrayList<LibraryEntry>> entriesBySource = getEntriesBySource(originalEntriesBySourceFile);
		
		ArrayList<LibraryEntry> corrected=new ArrayList<>();
		for (String source : entriesBySource.keySet()) {
			String sourceName=FilenameUtils.getBaseName(source);
			
			ArrayList<LibraryEntry> entries = entriesBySource.get(source);
			RetentionTimeFilter filter=getFilter(singleSource, entries, new File(libraryFile.getParent(), sourceName));
			Logger.logLine("Found "+filter.size()+" mutual data points with "+sourceName);

			for (LibraryEntry entry : entries) {
				LibraryEntry adjusted=entry.updateRetentionTime(filter.getXValue(entry.getScanStartTime()));
				corrected.add(adjusted);
			}
		}
		return corrected;
	}
	
	public static ArrayList<LibraryEntry> correctRTs(HashMap<String, ArrayList<LibraryEntry>> originalEntriesBySourceFile, File libraryFile) {
		HashMap<String, ArrayList<LibraryEntry>> entriesBySource = getEntriesBySource(originalEntriesBySourceFile);
		
		ArrayList<IndexedObject<String>> sourcesBySize=new ArrayList<>();
		for (Entry<String, ArrayList<LibraryEntry>> entry : entriesBySource.entrySet()) {
			sourcesBySize.add(new IndexedObject<String>(entry.getValue().size(), entry.getKey()));
		}
		Collections.sort(sourcesBySize);
		
		// get the biggest source of entries
		IndexedObject<String> startingDataset = sourcesBySize.remove(sourcesBySize.size()-1);
		ArrayList<LibraryEntry> corrected=new ArrayList<LibraryEntry>(entriesBySource.get(startingDataset.y));
		Logger.logLine("Starting alignment "+startingDataset+" with "+corrected.size()+" data points");
		
		// for each next source (in order)
		while (sourcesBySize.size()>0) {
			String source=sourcesBySize.remove(sourcesBySize.size()-1).y;
			String sourceName=FilenameUtils.getBaseName(source);
			
			ArrayList<LibraryEntry> entries = entriesBySource.get(source);
			RetentionTimeFilter filter=getFilter(corrected, entries, new File(libraryFile.getParent(), sourceName));

			for (LibraryEntry entry : entries) {
				LibraryEntry adjusted=entry.updateRetentionTime(filter.getXValue(entry.getScanStartTime()));
				corrected.add(adjusted);
			}
		}
		return corrected;
	}

	private static HashMap<String, ArrayList<LibraryEntry>> getEntriesBySource(
			HashMap<String, ArrayList<LibraryEntry>> originalEntriesBySourceFile) {
		HashMap<String, ArrayList<LibraryEntry>> entriesBySource=new HashMap<>();

		for (Entry<String, ArrayList<LibraryEntry>> originalEntries : originalEntriesBySourceFile.entrySet()) {
			for (LibraryEntry entry : originalEntries.getValue()) {
				ArrayList<LibraryEntry> list=entriesBySource.get(originalEntries.getKey());
				if (list==null) {
					list=new ArrayList<>();
					entriesBySource.put(originalEntries.getKey(), list);
				}
				list.add(entry);
			}
		}
		return entriesBySource;
	}

	/**
	 * longer and shorter lists are only for search speed considerations -- it's fine if they're reversed
	 * @param first preferably the longer list
	 * @param second preferably the shorter list
	 * @return
	 */
	protected static RetentionTimeFilter getFilter(ArrayList<LibraryEntry> first, ArrayList<LibraryEntry> second, File libraryFile) {
		HashMap<String, TFloatArrayList> rtsByPeptideModSeqFirst = getRTsBySeq(first);
		HashMap<String, TFloatArrayList> rtsByPeptideModSeqSecond = getRTsBySeq(second);
		String secondSource=second.get(0).getSource();
		
		return getFilter(rtsByPeptideModSeqFirst, rtsByPeptideModSeqSecond, libraryFile, secondSource);
	}

	private static RetentionTimeFilter getFilter(HashMap<String, TFloatArrayList> rtsByPeptideModSeqFirst, HashMap<String, TFloatArrayList> rtsByPeptideModSeqSecond, File libraryFile, String secondSource) {
		ArrayList<XYPoint> points=new ArrayList<>();
		for (Entry<String, TFloatArrayList> entry : rtsByPeptideModSeqFirst.entrySet()) {
			TFloatArrayList secondRTs=rtsByPeptideModSeqSecond.get(entry.getKey());
			if (secondRTs!=null) {
				points.add(new PeptideXYPoint(QuickMedian.median(entry.getValue().toArray()), QuickMedian.median(secondRTs.toArray()), false, entry.getKey()));
			}
		}
		
		Logger.logLine("Adding "+secondSource+" with "+points.size()+" matched data points");
		RetentionTimeFilter filter = RetentionTimeFilter.getFilter(points, "RT (seconds)", secondSource, TwoDimensionalKDE.HIGHER_RESOLUTION);
		filter.plot(points, Optional.ofNullable(libraryFile), "Global", secondSource);
		return filter;
	}

	protected static HashMap<String, TFloatArrayList> getRTsBySeq(ArrayList<LibraryEntry> list) {
		HashMap<String, TFloatArrayList> rtsByPeptideModSeq=new HashMap<>();
		for (LibraryEntry entry : list) {
			if (entry.isDecoy()) continue; // drop all decoys
			
			TFloatArrayList rts=rtsByPeptideModSeq.get(entry.getPeptideModSeq());
			if (rts==null) {
				rts=new TFloatArrayList();
				rtsByPeptideModSeq.put(entry.getPeptideModSeq(), rts);
			}
			rts.add(entry.getScanStartTime());
		}
		return rtsByPeptideModSeq;
	}
	
	public static ArrayList<LibraryEntry> filterIons(ArrayList<LibraryEntry> entries, float minimumIntensityFraction) {
		ArrayList<LibraryEntry> processed=new ArrayList<>();
		for (LibraryEntry entry : entries) {
			double[] massArray=entry.getMassArray();
			float[] intensityArray=entry.getIntensityArray();
			float[] correlationArray=entry.getCorrelationArray();
			boolean[] quantifiedIonsArray=entry.getQuantifiedIonsArray();
			ArrayList<PeakChromatogram> peaks=new ArrayList<>();
			int numPeaks=Math.min(massArray.length, correlationArray.length);
			float minimumIntensity=minimumIntensityFraction*General.max(intensityArray);
			for (int i=0; i<numPeaks; i++) {
				if (intensityArray[i]>minimumIntensity) {
					peaks.add(new PeakChromatogram(massArray[i], intensityArray[i], correlationArray[i], quantifiedIonsArray[i]));
				}
			}
			Collections.sort(peaks);
			Quadruplet<double[], float[], float[], boolean[]> arrays=PeakChromatogram.toChromatogramArrays(peaks);
			
			massArray=arrays.x;
			intensityArray=arrays.y;
			correlationArray=arrays.z;
			quantifiedIonsArray=arrays.w;
			processed.add(entry.updateMS2(massArray, intensityArray, correlationArray, quantifiedIonsArray));
		}
		return processed;
	}
}
