package edu.washington.gs.maccoss.encyclopedia.algorithms.curve;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.curve.DilutionCurveFitter.AlignmentWithAnchors;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakIntensityComparator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class TargetedScheduler {
	
	public static void main(String[] args) throws Exception {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		final File outputDirectory=new File("/Users/searle.30/Documents/CCIC/maisam/112922_pancreatitis_targeted_dataset/urine_29nov2022/wide_scheduled/");
		File libraryFile=new File("/Users/searle.30/Documents/CCIC/maisam/112922_pancreatitis_targeted_dataset/urine_29nov2022/112922_urine_pool_GPFDIA_2.5ul_combined.dia.elib");
		File rtAlignFile=new File("/Users/searle.30/Documents/CCIC/maisam/112922_pancreatitis_targeted_dataset/urine_29nov2022/112922_urine_pool_DIA_05.dia.elib");
		File targetFastaFile=new File("/Users/searle.30/Documents/CCIC/maisam/112922_pancreatitis_targeted_dataset/urine_29nov2022/urine_target_proteins_for_AP_and_CP.fasta");
		HashSet<String> keyAccessionNumbers=new HashSet<>();
		keyAccessionNumbers.add("sp|P04746|AMYP_HUMAN");
		keyAccessionNumbers.add("sp|P08217|CEL2A_HUMAN");
		keyAccessionNumbers.add("sp|P02741|CRP_HUMAN");
		keyAccessionNumbers.add("sp|Q8WVV4|POF1B_HUMAN");
		keyAccessionNumbers.add("sp|P05451|REG1A_HUMAN");
		keyAccessionNumbers.add("sp|P04745|AMY1_HUMAN");
		keyAccessionNumbers.add("sp|P02765|FETUA_HUMAN");
		keyAccessionNumbers.add("sp|A0A075B6J9|LV218_HUMAN");
		keyAccessionNumbers.add("sp|P02647|APOA1_HUMAN");
		keyAccessionNumbers.add("sp|P19961|AMY2B_HUMAN");
		keyAccessionNumbers.add("sp|P02671|FIBA_HUMAN");
		keyAccessionNumbers.add("sp|P02675|FIBB_HUMAN");
		keyAccessionNumbers.add("sp|P02679|FIBG_HUMAN");
		keyAccessionNumbers.add("sp|P02652|APOA2_HUMAN");
		keyAccessionNumbers.add("sp|P35030|TRY3_HUMAN");
		keyAccessionNumbers.add("sp|P08246|ELNE_HUMAN");
		keyAccessionNumbers.add("sp|P07478|TRY2_HUMAN");
		
		AbstractDilutionCurveFittingParameters fittingParams=new Targeted10HzParametersWithWideWindows();
		generateAssay(params, outputDirectory, libraryFile, rtAlignFile, targetFastaFile, keyAccessionNumbers, fittingParams);
	}
	
	public static void main2(String[] args) throws Exception {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		final File outputDirectory=new File("/Users/searleb/Documents/OSU/projects/yi/051622/scheduled/");
		File libraryFile=new File("/Users/searleb/Documents/OSU/projects/yi/051622/051622_Mouse_Tcell_pool_clib.elib");
		File rtAlignFile=new File("/Users/searleb/Documents/OSU/projects/yi/051622/051622_Mouse_Tcell_pool_DIA_07.mzML.elib");
		File targetFastaFile=new File("/Users/searleb/Documents/OSU/projects/yi/051622/immuno_oncology_targets.fasta");
		HashSet<String> keyAccessionNumbers=new HashSet<>();

		AbstractDilutionCurveFittingParameters fittingParams=new Targeted10HzParameters();
		generateAssay(params, outputDirectory, libraryFile, rtAlignFile, targetFastaFile, keyAccessionNumbers, fittingParams);
	}


	private final static PeakIntensityComparator intensityComparator=new PeakIntensityComparator();

	public static void generateAssay(SearchParameters params, final File outputDirectory, File libraryFile, File rtAlignFile, File targetFastaFile, HashSet<String> keyAccessionNumbers, AbstractDilutionCurveFittingParameters fittingParams) throws IOException, SQLException, DataFormatException, FileNotFoundException, UnsupportedEncodingException {
		final File exportLibraryFile=new File(outputDirectory, "target_library.dlib");
		final File libraryAlignmentFile=new File(outputDirectory, "library_rt_alignment.pdf");

		outputDirectory.mkdirs();

		final HashMap<String, LibraryEntry> libraryEntryByPeptideModSeq=DilutionCurveFitter.getLibraryData(params, libraryFile);
		
		final TObjectFloatHashMap<String> knownRTInSecs=new TObjectFloatHashMap<String>();
		ArrayList<XYPoint> rts=new ArrayList<XYPoint>();
		for (Entry<String, LibraryEntry> entry : DilutionCurveFitter.getLibraryData(params, rtAlignFile).entrySet()) {
			LibraryEntry idealEntry=libraryEntryByPeptideModSeq.get(entry.getKey());
			float alignmentRT = entry.getValue().getScanStartTime();
			knownRTInSecs.put(entry.getKey(), alignmentRT);
			if (idealEntry!=null) {
				XYPoint xy=new XYPoint(idealEntry.getScanStartTime()/60f, alignmentRT/60f);
				rts.add(xy);
			}
		}
		
		RetentionTimeFilter rtAlignmentFilter=RetentionTimeFilter.getFilter(rts, "Library Retention Time (min)", "Alignment Retention Time (min)");
		rtAlignmentFilter.plot(rts, Optional.of(libraryAlignmentFile));
		final AlignmentWithAnchors rtAlignment=new AlignmentWithAnchors(rtAlignmentFilter, knownRTInSecs);

		float minRTInSec=Float.MAX_VALUE;
		float maxRTInSec=-Float.MAX_VALUE;
		for (LibraryEntry entry : libraryEntryByPeptideModSeq.values()) {
			float rtInSec = rtAlignment.getAlignedRTInSec(entry, true);
			if (rtInSec>maxRTInSec) maxRTInSec=rtInSec;
			if (rtInSec<minRTInSec) minRTInSec=rtInSec;
		}
		Range rtInSecRange=new Range(minRTInSec, maxRTInSec);
		//rtInSecRange=new Range(12*60f, 95*60f);
		
		ArrayList<FastaEntryInterface> targetProteins=FastaReader.readFasta(targetFastaFile, params);
		HashSet<String> targetAccessionNumbers=new HashSet<>();
		for (FastaEntryInterface entry : targetProteins) {
			targetAccessionNumbers.add(entry.getAccession());
		}
		targetAccessionNumbers.removeAll(keyAccessionNumbers);
		
		ArrayList<ScoredEntry> keyTargetEntries=new ArrayList<ScoredEntry>();
		ArrayList<ScoredEntry> optionalTargetEntries=new ArrayList<ScoredEntry>();
		for (LibraryEntry entry : libraryEntryByPeptideModSeq.values()) {

			ArrayList<Peak> peaks=entry.getPeaks(TransitionRefiner.identificationCorrelationThreshold);
			Collections.sort(peaks, intensityComparator);
			Collections.reverse(peaks);
			
			float top3Ions=0;
			float top6Ions=0;
			int peakNum=0;
			for (Peak peak : peaks) {
				peakNum++;
				if (peakNum==3) {
					top3Ions=peak.intensity;
				}
				if (peakNum==6) {
					top6Ions=peak.intensity;
				}
			}
			
			for (String accession : entry.getAccessions()) {
				
				if (keyAccessionNumbers.contains(accession)) {
					keyTargetEntries.add(new ScoredEntry(top3Ions, entry));
					break;
				} else if (targetAccessionNumbers.contains(accession)) {
					optionalTargetEntries.add(new ScoredEntry(top3Ions, entry));
					break;
				}
			}
		}
		Collections.sort(keyTargetEntries);
		Collections.reverse(keyTargetEntries);
		Collections.sort(optionalTargetEntries);
		Collections.reverse(optionalTargetEntries);
		
		ArrayList<LibraryEntry> potentialTargetEntries=new ArrayList<>();
		for (ScoredEntry e : keyTargetEntries) {
			potentialTargetEntries.add(e.y);
		}
		for (ScoredEntry e : optionalTargetEntries) {
			potentialTargetEntries.add(e.y);
		}

		ArrayList<LibraryEntry> targetEntries = scheduleAssay(outputDirectory, 
				fittingParams, libraryEntryByPeptideModSeq, rtAlignment, rtInSecRange, potentialTargetEntries);
		
		HashSet<String> allAccessionNumbers=new HashSet<String>();
		allAccessionNumbers.addAll(keyAccessionNumbers);
		allAccessionNumbers.addAll(targetAccessionNumbers);
		for (LibraryEntry entry : targetEntries) {
			HashSet<String> accessions=entry.getAccessions();
			for (String string : accessions) {
				allAccessionNumbers.remove(string);
			}
		}

		DilutionCurveFitter.writeLibraryEntries(params, exportLibraryFile, targetEntries);
		
		Logger.logLine("Failed to schedule "+allAccessionNumbers.size()+" proteins:");
		for (String string : allAccessionNumbers) {
			Logger.logLine("\t"+string);
		}
	}



	private static ArrayList<LibraryEntry> scheduleAssay(final File outputDirectory, AbstractDilutionCurveFittingParameters fittingParams,
			final HashMap<String, LibraryEntry> libraryEntryByPeptideModSeq, final AlignmentWithAnchors rtAlignment, 
			Range rtInSecRange, final ArrayList<LibraryEntry> potentialTargetEntries)
			throws FileNotFoundException, UnsupportedEncodingException {

		ArrayList<LibraryEntry> targetEntries=new ArrayList<LibraryEntry>();
		boolean hitMaxDensity=false;
		float[] assayRT=new float[Math.round(rtInSecRange.getStop()+fittingParams.getWindowInMin()*60f)]; // N+W minutes in second increments
		for (int i = 0; i < assayRT.length; i++) {
			assayRT[i]=i/60f;
		}
		float[] assayDensity=new float[assayRT.length];
		final PrintWriter assayWriter=new PrintWriter(new File(outputDirectory, "assay.csv"), "UTF-8");
		assayWriter.println("Compound,Formula,Adduct,m/z,z,RT Time (min),Window (min)");
		
		int count=0;
		HashMap<String, ArrayList<LibraryEntry>> targetPeptidesByProtein=new HashMap<String, ArrayList<LibraryEntry>>();
		for (LibraryEntry entry : potentialTargetEntries) {
			String accessionsKey=PSMData.accessionsToString(entry.getAccessions());
			if (entry.getAccessions().size()>1) continue;
			ArrayList<LibraryEntry> list=targetPeptidesByProtein.get(accessionsKey);
			if (list==null) {
				list=new ArrayList<LibraryEntry>();
				targetPeptidesByProtein.put(accessionsKey, list);
			}

			boolean keep=true;
			if (count<fittingParams.getTargetTotalNumberOfPeptides()) {
				if (list.size()<fittingParams.getMaxNumberPeptidesPerProtein()) {

					float rtInSec = rtAlignment.getAlignedRTInSec(entry, true);
					float[] testDensity=DilutionCurveFitter.incrementDensity(rtInSec, fittingParams.getWindowInMin(rtInSec), assayDensity);
					for (int i = 0; i < testDensity.length; i++) {
						if (testDensity[i]>fittingParams.getAssayMaxDensity()) {
							keep=false;
							
							if (!hitMaxDensity) {
								hitMaxDensity=true;
								Logger.logLine("First hit of max density at LOQ: "+entry.getScore());
							}
							break;
						}
					}
					
					if (keep) {
						assayDensity=testDensity; // update density
						count++;
						Logger.logLine("Adding peptide ("+count+") to assay: ("+(rtAlignment.isKnown(entry)?" - ":"!!!")+") "+entry.getPeptideModSeq()+" --> PEP: "+entry.getScore()+" from "+accessionsKey);
						list.add(entry);
					}
				}
			}
		}

		count=0;
		int numSingletons=0;
		ArrayList<String> keys=new ArrayList<String>(targetPeptidesByProtein.keySet());
		Collections.sort(keys);
		for (String key : keys) {
			ArrayList<LibraryEntry> list=targetPeptidesByProtein.get(key);
			for (LibraryEntry entry : list) {
				float rtInSec = rtAlignment.getAlignedRTInSec(entry, false);
				DilutionCurveFitter.addPeptideToAssay(assayWriter, entry, rtInSec, fittingParams.getWindowInMin(rtInSec));
				targetEntries.add(entry.updateRetentionTime(rtInSec));
				count++;
			}
			if (list.size()==1) {
				numSingletons++;
			} if (list.size()==0) {
				targetPeptidesByProtein.remove(key);
			}
		}

		
		assayWriter.flush();
		assayWriter.close();
		Logger.logLine("Finished writing assay for "+targetPeptidesByProtein.size()+" proteins using "+count+" total peptides ("+numSingletons+" single peptide targets)");
		DilutionCurveFitter.writeSchedulingGraph(outputDirectory, assayRT, assayDensity);
				
		return targetEntries;
	}
	public static class ScoredEntry extends Pair<Float, LibraryEntry> implements Comparable<ScoredEntry> {

		public ScoredEntry(float x, LibraryEntry y) {
			super(x, y);
		}

		@Override
		public int compareTo(ScoredEntry o) {
			if (o==null) return 1;
			int c=Float.compare(x, o.x);
			if (c!=0) return c;
			
			return y.compareTo(o.y);
		}
		
		public float getScore() {
			return x;
		}
	}

}
