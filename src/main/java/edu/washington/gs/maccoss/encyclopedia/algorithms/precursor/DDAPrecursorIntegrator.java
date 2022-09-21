package edu.washington.gs.maccoss.encyclopedia.algorithms.precursor;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMScoredSpectrum;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.RetentionTimeComparator;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;

public class DDAPrecursorIntegrator {
	private static final Range targetMzRange=new Range(-Float.MAX_VALUE, Float.MAX_VALUE);
	public static ArrayList<IntegratedLibraryEntry> integrateSearch(ProgressIndicator progress, HashMap<String, PSMData> targetPSMs, StripeFileInterface stripeFile, final SearchParameters parameters) {
		
		HashMap<String, ArrayList<PSMScoredSpectrum>> psmsBySequence=new HashMap<>();
		
		float rtTotalWindowSize = parameters.getExpectedPeakWidth()*8; // pretty wide, better to consider more in case we're at a shoulder 
		ArrayList<PSMData> psms=new ArrayList<PSMData>(targetPSMs.values());
		Collections.sort(psms, RetentionTimeComparator.comparator);
		float minRT=psms.get(0).getRetentionTimeInSec()-rtTotalWindowSize/2f;
		float maxRT=psms.get(psms.size()-1).getRetentionTimeInSec()+rtTotalWindowSize/2;
		float increment=(maxRT-minRT)/10f;
		

		try {
			float currentRTStart=minRT;
			float currentRTStop=minRT+increment;
			ArrayList<FragmentScan> scanBlock=stripeFile.getStripes(targetMzRange, currentRTStart, currentRTStop, false);
			for (PSMData psm : psms) {
				if (psm.getRetentionTimeInSec()>currentRTStop) {
					Logger.logLine("Reading MS/MS from "+stripeFile.getOriginalFileName()+", from "+(currentRTStart/60)+" to "+(currentRTStop/60)+" mins");
					currentRTStart=currentRTStop;
					currentRTStop=currentRTStart+increment;
					scanBlock=stripeFile.getStripes(targetMzRange, currentRTStart, currentRTStop, false);
				}
				
				ArrayList<PSMScoredSpectrum> list=psmsBySequence.get(psm.getPeptideModSeq());
				if (list==null) {
					list=new ArrayList<>();
					psmsBySequence.put(psm.getPeptideModSeq(), list);
				}
				
				FragmentScan matchingMSMS=null;
				for (FragmentScan msms : scanBlock) {
					if (msms.getSpectrumIndex()==psm.getSpectrumIndex()) {
						matchingMSMS=msms;
						break;
					}
				}
				
				if (matchingMSMS!=null) {
					list.add(new PSMScoredSpectrum(psm, matchingMSMS));
				} else {
					// create a "fake" scan at this RT if one wasn't collected so that we can quantify the precursor
					FragmentScan emptyScan=new FragmentScan("N/A", "N/A", 0, psm.getRetentionTimeInSec(), 0, 0.0f, psm.getPrecursorMZ()-1, psm.getPrecursorMZ()+1,
							new double[0], new float[0]);
					list.add(new PSMScoredSpectrum(psm, emptyScan));
				}
				
			}
			
			final PrecursorIonMap map=new PrecursorIonMap(parameters.getPrecursorTolerance(), rtTotalWindowSize);
			for (Entry<String, ArrayList<PSMScoredSpectrum>> entry : psmsBySequence.entrySet()) {
				map.addPeptide(entry.getKey(), stripeFile.getOriginalFileName(), entry.getValue());
			}
			map.finalizeMap();
			
			SubProgressIndicator sub1=new SubProgressIndicator(progress, 0.4f);
			Logger.logLine("Reading precursors from "+stripeFile.getOriginalFileName());
			ArrayList<PrecursorScan> precursors = stripeFile.getPrecursors(minRT, maxRT);
			int count=0;
			for (PrecursorScan scan : precursors) {
				map.processSpectrum(scan);
				count++;
				if (count%1000==0) {
					String message = "Processed "+(scan.getScanStartTime()/60f)+" mins of precursors";
					sub1.update(message, count/(float)precursors.size());
					Logger.logLine(message);
				}
			}

			SubProgressIndicator sub2=new SubProgressIndicator(progress, 0.5f);
			Logger.logLine("Integrating peaks from "+stripeFile.getOriginalFileName());

			ArrayList<IntegratedPeptide> peptides=map.getPeptides();
			ArrayList<IntegratedLibraryEntry> entries=new ArrayList<>();
			for (IntegratedPeptide peptide : peptides) {
				entries.addAll(peptide.integrate(parameters));
			}
			return entries;
		
		} catch (IOException | SQLException | DataFormatException e) {
			Logger.errorException(e);
			e.printStackTrace();
			throw new EncyclopediaException("Failed to parse DDA data for quantification", e);
		}
	}
}
