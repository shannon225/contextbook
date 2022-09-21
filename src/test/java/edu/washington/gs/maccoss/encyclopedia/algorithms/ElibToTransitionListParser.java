package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakIntensityComparator;

public class ElibToTransitionListParser {
	private final static PeakIntensityComparator intensityComparator=new PeakIntensityComparator();
	
	public static void main(String[] args) throws Exception {
		File elib=new File("/Users/searle.30/Documents/prosit/hela_specific_clib.z2_nce33.elib");
		LibraryFile library=new LibraryFile();
		library.openFile(elib);

		HashMap<String, ArrayList<LibraryEntry>> proteinMap=new HashMap<>();
		ArrayList<LibraryEntry> entries = library.getAllEntries(false, new AminoAcidConstants());
		for (LibraryEntry entry : entries) {
			ArrayList<Peak> peaks=entry.getPeaks(TransitionRefiner.identificationCorrelationThreshold);
			Collections.sort(peaks, intensityComparator);
			Collections.reverse(peaks);
			
			float top3Ions=0;
			int peakNum=0;
			for (Peak peak : peaks) {
				peakNum++;
				if (peakNum==3) {
					top3Ions=peak.intensity;
				}
			}
			// require at least 3 ions!
			if (top3Ions==0) continue;
			
			String proteins=PSMData.accessionsToString(entry.getAccessions());
			
			ArrayList<LibraryEntry> list=proteinMap.get(proteins);
			if (list==null) {
				list=new ArrayList<LibraryEntry>();
				proteinMap.put(proteins, list);
			}
			list.add(entry);
		}

		System.out.println("ProteinIndex\tNumPeptides\tPeptideModSeq\tPeptideSeq\tPrecursorCharge\tRetentionTimeInSec\t3rdIonIntensity\t6thIonIntensity\tTIC");
		for (Entry<String, ArrayList<LibraryEntry>> mapping : proteinMap.entrySet()) {
			ArrayList<LibraryEntry> list=mapping.getValue();
			if (list.size()<=1) continue;
			
			String proteins=mapping.getKey();
			Collections.sort(list);
			
			for (LibraryEntry entry : entries) {
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
				System.out.println(proteins+"\t"+list.size()+"\t"+entry.getPeptideModSeq()+"\t"+entry.getPeptideSeq()+"\t"+entry.getPrecursorCharge()+"\t"+entry.getRetentionTime()+"\t"+top3Ions+"\t"+top6Ions+"\t"+entry.getTIC());
			}
		}
	}

}
