package edu.washington.gs.maccoss.encyclopedia.algorithms.curve;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.curve.DilutionCurveFitter.AlignmentWithAnchors;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class TargetedScheduler {
	
	public static void main(String[] args) throws Exception {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		final File outputDirectory=new File("/Users/searleb/Documents/OSU/projects/yi/051622/scheduled/");
		outputDirectory.mkdirs();
		File libraryFile=new File("/Users/searleb/Documents/OSU/projects/yi/051622/051622_Mouse_Tcell_pool_clib.elib");
		File rtAlignFile=new File("/Users/searleb/Documents/OSU/projects/yi/051622/051622_Mouse_Tcell_pool_DIA_07.mzML.elib");
		File targetFastaFile=new File("/Users/searleb/Documents/OSU/projects/yi/051622/immuno_oncology_targets.fasta");
		final File exportLibraryFile=new File(outputDirectory, "target_library.dlib");
		
		final File libraryAlignmentFile=new File(outputDirectory, "library_rt_alignment.pdf");
		
		AbstractDilutionCurveFittingParameters fittingParams=new Targeted10HzParameters();

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
			float rtInSec = rtAlignment.getAlignedRTInSec(entry);
			if (rtInSec>maxRTInSec) maxRTInSec=rtInSec;
			if (rtInSec<minRTInSec) minRTInSec=rtInSec;
		}
		Range rtInSecRange=new Range(minRTInSec, maxRTInSec);
		rtInSecRange=new Range(12*60f, 95*60f);
		
		ArrayList<FastaEntryInterface> targetProteins=FastaReader.readFasta(targetFastaFile, params);
		HashSet<String> targetAccessionNumbers=new HashSet<>();
		for (FastaEntryInterface entry : targetProteins) {
			targetAccessionNumbers.add(entry.getAccession());
		}
		ArrayList<LibraryEntry> potentialTargetEntries=new ArrayList<LibraryEntry>();
		for (LibraryEntry entry : libraryEntryByPeptideModSeq.values()) {
			for (String accession : entry.getAccessions()) {
				if (targetAccessionNumbers.contains(accession)) {
					potentialTargetEntries.add(entry);
					break;
				}
			}
		}
		Collections.sort(potentialTargetEntries, new Comparator<LibraryEntry>() {
			@Override
			public int compare(LibraryEntry o1, LibraryEntry o2) {
				if (o1==null&&o2==null) return 0;
				if (o1==null) return 1;
				if (o2==null) return -1;
				if(o1.getScore()>o2.getScore()) {
					return 1; // low score is better
				}
				if (o1.getScore()<o2.getScore()) {
					return -1;
				}
				return o1.compareTo(o2);
			}
		});

		ArrayList<LibraryEntry> targetEntries = scheduleAssay(outputDirectory, 
				fittingParams, libraryEntryByPeptideModSeq, rtAlignment, rtInSecRange, potentialTargetEntries);

		DilutionCurveFitter.writeLibraryEntries(params, exportLibraryFile, targetEntries);
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

					float rtInSec = rtAlignment.getAlignedRTInSec(entry);
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
						Logger.logLine("Adding peptide ("+count+") to assay: "+entry.getPeptideModSeq()+" --> PEP: "+entry.getScore()+" from "+accessionsKey);
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
				float rtInSec = rtAlignment.getAlignedRTInSec(entry);
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
}
