package edu.washington.gs.maccoss.encyclopedia.algorithms.precursor;

import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.HasRetentionTime;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMScoredSpectrum;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeakTrace;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakWithTime;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PrecursorIon;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TByteDoubleHashMap;
import gnu.trove.procedure.TByteDoubleProcedure;

public class IntegratedPeptide implements HasRetentionTime {
	public static final byte[] integratedIsotopes=new byte[] {-1, 0, 1, 2};
	private final String peptideModSeq;
	private final String sourceFile;
	private final ArrayList<PSMScoredSpectrum> psms;
	private final float medianRT;
	private final PrecursorIon[] ions;
	private final ArrayList<PeakWithTime>[] peaks; // mutable!
	
	@SuppressWarnings("unchecked")
	public IntegratedPeptide(String peptideModSeq, String sourceFile, ArrayList<PSMScoredSpectrum> psms) {
		this.peptideModSeq = peptideModSeq;
		this.sourceFile=sourceFile;
		this.psms = psms;
		TFloatArrayList rts=new TFloatArrayList();
		TByteDoubleHashMap precursors=new TByteDoubleHashMap();
		for (PSMScoredSpectrum psm : psms) {
			rts.add(psm.getRetentionTimeInSec());
			precursors.put(psm.getCharge(), psm.getMZ());
		}
		medianRT=QuickMedian.median(rts.toArray());

		ArrayList<PrecursorIon> ionList=new ArrayList<>();
		precursors.forEachEntry(new TByteDoubleProcedure() {
			@Override
			public boolean execute(byte charge, double mz) {
				for (int i = 0; i < integratedIsotopes.length; i++) {
					ionList.add(new PrecursorIon(peptideModSeq, mz+integratedIsotopes[i]*MassConstants.neutronMass/charge, charge, integratedIsotopes[i]));
				}
				return true;
			}
		});
		
		ions=ionList.toArray(new PrecursorIon[ionList.size()]);
		Arrays.sort(ions);	
		
		peaks=new ArrayList[ions.length];
		for (int i = 0; i < ions.length; i++) {
			peaks[i]=new ArrayList<>();
		}
	}
	
	@Override
	public float getRetentionTimeInSec() {
		return medianRT;
	}
	
	public ArrayList<PSMScoredSpectrum> getPsms() {
		return psms;
	}
	
	public String getKey() {
		return peptideModSeq;
	}
	
	public PrecursorIon[] getIons() {
		return ions;
	}
	
	public ArrayList<PeakWithTime>[] getPeaks() {
		return peaks;
	}

	public ArrayList<PeakTrace<PrecursorIon>> getCenteredTraces() {
		return getTraces(true, true);
	}
	private ArrayList<PeakTrace<PrecursorIon>> getTraces(boolean normalize, boolean center) {
		ArrayList<PeakTrace<PrecursorIon>> traces=new ArrayList<>();
		for (int i = 0; i < ions.length; i++) {
			ArrayList<PeakWithTime> peakList = peaks[i];
			PeakTrace<PrecursorIon> trace = getTrace(ions[i], peakList, normalize, center);
			traces.add(trace);
		}
		
		return traces;
	}

	private static PeakTrace<PrecursorIon> getTrace(PrecursorIon ion, ArrayList<PeakWithTime> peakList, boolean normalize, boolean center) {
		float maxIntensity=-1;
		float maxRT=0;
		for (PeakWithTime peak : peakList) {
			if (maxIntensity<peak.getIntensity()) {
				maxIntensity=peak.getIntensity();
				maxRT=peak.getRtInSec();
			}
		}
		if (!center) maxRT=0;
		if (!normalize) maxIntensity=1.0f;
		
		TFloatArrayList x=new TFloatArrayList();
		TFloatArrayList y=new TFloatArrayList();

		for (PeakWithTime peak : peakList) {
			x.add(peak.getRtInSec()-maxRT);
			y.add(peak.getIntensity()/maxIntensity);
		}
		PeakTrace<PrecursorIon> trace = new PeakTrace<PrecursorIon>(ion, maxRT, x.toArray(), y.toArray());
		return trace;
	}
	
	public ArrayList<IntegratedLibraryEntry> integrate(SearchParameters params) {
		ArrayList<PeakTrace<PrecursorIon>> traces=getTraces(false, false);
		TFloatArrayList rts=new TFloatArrayList();
		for (PSMScoredSpectrum pep : psms) {
			rts.add(pep.getRetentionTimeInSec());
		}
		
		float[] isotopicDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(psms.get(0).getPeptideData().getPeptideModSeq(), params.getAAConstants());
		ArrayList<TransitionRefinementData> trd=PrecursorIntegrator.integratePeptide(peptideModSeq, isotopicDistribution, rts.toArray(), traces, params);
		
		ArrayList<IntegratedLibraryEntry> entries=new ArrayList<>();
		for (TransitionRefinementData data : trd) {
			PSMScoredSpectrum bestPSM=null;
			int copies=0;
			for (PSMScoredSpectrum pep : psms) {
				if (pep.getCharge()==data.getPrecursorCharge()) {
					copies++;
					if (bestPSM==null||bestPSM.getPeptideData().getSortingScore()<pep.getPeptideData().getSortingScore()) {
						bestPSM=pep;
					}
				}
			}

			if (bestPSM==null) {
				Logger.errorLine("Missing PSM for "+data.getPeptideModSeq()+" +"+data.getPrecursorCharge()+"H");
			} else {
				PSMData psm=bestPSM.getPeptideData();
				FragmentScan msms=bestPSM.getMsms();
				float[] fakeCorrelationArray=new float[msms.getMassArray().length];
				Arrays.fill(fakeCorrelationArray, 1.0f);
				IntegratedLibraryEntry entry=new IntegratedLibraryEntry(sourceFile, psm.getAccessions(), psm.getSpectrumIndex(), psm.getPrecursorMZ(), psm.getPrecursorCharge(), 
						psm.getPeptideModSeq(), copies, msms.getScanStartTime(), psm.getScore(), msms.getMassArray(), msms.getIntensityArray(), fakeCorrelationArray, data);
				entries.add(entry);
			}
		}
		
		return entries;
	}
	
	public String toString() {
		ArrayList<PeakTrace<PrecursorIon>> traces=getTraces(false, false);

		StringBuilder sb=new StringBuilder("ArrayList<PeakTrace<Ion>> traces=new ArrayList<>();");
		int index=0;
		for (PeakTrace<PrecursorIon> trace : traces) {
			index++;
			PrecursorIon ion=trace.getIon();
			sb.append("PrecursorIon ion"+index+"=new PrecursorIon(\""+ion.getName()+"\", "+ion.getMass()+", (byte)"+ion.getCharge()+", (byte)"+ion.getIsotope()+");");
			sb.append("float[] rt"+index+"=new float[] {"+General.toString(trace.getRt(), "f,")+"f};");
			sb.append("float[] intensity"+index+"=new float[] {"+General.toString(trace.getIntensity(), "f,")+"f};");
			sb.append("PeakTrace<PrecursorIon> trace"+index+"=new PeakTrace<Ion>(ion"+index+", "+trace.getRetentionTimeInSec()+"f, rt"+index+", intensity"+index+");");
			sb.append("traces.add(trace"+index+");");
		}
		return sb.toString();
		
	}
}
