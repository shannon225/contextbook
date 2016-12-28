package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Optional;
import java.util.concurrent.ConcurrentLinkedQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.AmbiguousPeptideModSeq;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYZPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SkylineSGFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ThreadableTask;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class PeptideQuantExtractorTask extends ThreadableTask<Nothing> {
	private final Optional<PhosphoLocalizer> localizer;
	private final String filename;
	private final ArrayList<Stripe> stripes;
	private final boolean limitToQuantifiable;

	private final PSMScorer scorer;
	private final SearchParameters params;

	private final PSMData psmdata;
	private final ConcurrentLinkedQueue<IntegratedLibraryEntry> savedEntries; // CAN BE NULL

	public PeptideQuantExtractorTask(String filename, PSMData psmdata, Optional<PhosphoLocalizer> localizer, ArrayList<Stripe> stripes, SearchParameters parameters, boolean limitToQuantifiable) {
		this.filename=filename;
		this.psmdata=psmdata;
		this.localizer=localizer;
		this.stripes=stripes;

		scorer=new DotProduct(parameters);
		params=parameters;
		this.savedEntries=null;
		
		this.limitToQuantifiable=limitToQuantifiable; //library.isPresent();
	}

	public PeptideQuantExtractorTask(String filename, PSMData psmdata, Optional<PhosphoLocalizer> localizer, ArrayList<Stripe> stripes, SearchParameters parameters, ConcurrentLinkedQueue<IntegratedLibraryEntry> savedEntries, boolean limitToQuantifiable) {
		this.filename=filename;
		this.psmdata=psmdata;
		this.localizer=localizer;
		this.stripes=stripes;

		scorer=new DotProduct(parameters);
		params=parameters;
		this.savedEntries=savedEntries;
		
		this.limitToQuantifiable=limitToQuantifiable; //library.isPresent();
	}
	
	public ArrayList<Stripe> getScanSubset(float minRT, float maxRT) {
		ArrayList<Stripe> subset=new ArrayList<Stripe>();
		for (Stripe scan : stripes) {
			if (scan.getScanStartTime()>=minRT&&scan.getScanStartTime()<=maxRT) {
				subset.add(scan);
			}
		}
		return subset;
	}
	
	@Override
	public String getTaskName() {
		return psmdata.getPeptideModSeq();
	}

	@Override
	protected Nothing process() {
		Optional<TransitionRefinementData> spectrum=extractSpectrum(psmdata.getAccessions(), psmdata.getPrecursorCharge(), psmdata.getPeptideModSeq(), psmdata.getRetentionTime(), psmdata.getDuration(), limitToQuantifiable);
		Optional<HashMap<String, TransitionRefinementData>> phosphoData=Optional.empty();
		if (canRunLocalization()) {
			Optional<PhosphoLocalizationData> localizationData=runLocalization();
			if (localizationData.isPresent()) {
				phosphoData=Optional.ofNullable(localizationData.get().getPassingForms());
			}
		}
		if (spectrum.isPresent()) {
			// FIXME need to not add duplicates!!!! for now just run SQL:
			// delete from entries where RowId not in (SELECT MIN(RowId) FROM entries GROUP BY PeptideModSeq, PrecursorCharge)
			TransitionRefinementData data=spectrum.get();
			data.setModificationQuantData(phosphoData);
			if (canRunLocalization()) {
				if (!phosphoData.isPresent()) {
					// no need to localize since there's only one form, so annotate this directly on the data object
					int numberOfMods=PeptideUtils.getNumberOfMods(psmdata.getPeptideModSeq(), AmbiguousPeptideModSeq.NOMINAL_MASS);
					if (numberOfMods>0) {
						FragmentIon[] localizingIons=new FragmentIon[0];
						float localizationScore=1000.0f;
						data.setModificationLocalizationData(Optional.of(new ModificationLocalizationData(AmbiguousPeptideModSeq.getUnambigous(psmdata.getPeptideModSeq(), params.getAAConstants()), data.getApexRT(), localizationScore, numberOfMods, true, localizingIons)));
					}
				} else if (phosphoData.get().size()==0) {
					// no confident localizations, so annotate this directly on the data object
					int numberOfMods=PeptideUtils.getNumberOfMods(psmdata.getPeptideModSeq(), AmbiguousPeptideModSeq.NOMINAL_MASS);
					if (numberOfMods>0) {
						FragmentIon[] localizingIons=new FragmentIon[0];
						float localizationScore=0.0f;
						data.setModificationLocalizationData(Optional.of(new ModificationLocalizationData(AmbiguousPeptideModSeq.getUnambigous(psmdata.getPeptideModSeq(), params.getAAConstants()), data.getApexRT(), localizationScore, numberOfMods, false, localizingIons)));
					}
				}
			}
					
			double[] fragmentMassArray=FragmentIon.getMasses(data.getFragmentMassArray());
			IntegratedLibraryEntry entry=new IntegratedLibraryEntry(filename, psmdata.getAccessions(), psmdata.getSpectrumIndex(), psmdata.getPrecursorMZ(), psmdata.getPrecursorCharge(), psmdata.getPeptideModSeq(), 1, psmdata.getRetentionTime(), psmdata.getScore(), fragmentMassArray, data.getIntegrationArray(), data);
			if (limitToQuantifiable) {
				if (entry.getIonCount()<4||entry.getTIC()<1.0f) {
					return Nothing.NOTHING;
				}
			}
			if (savedEntries!=null) {
				savedEntries.add(entry);
			}
		}
		return Nothing.NOTHING;
	}

	private boolean canRunLocalization() {
		return params.isRunPhosphoLocalization()&&localizer.isPresent();
	}

	public Optional<PhosphoLocalizationData> runLocalization() {
		return localizer.get().runDIAPhosphoLocalization(psmdata, stripes);
	}

	private Optional<TransitionRefinementData> extractSpectrum(HashSet<String> accessions, byte precursorCharge, String peptideModSeq, float retentionTime, float duration, boolean limitToQuantifiable) {
		FragmentationModel model=new FragmentationModel(peptideModSeq, params.getAAConstants());
		AnnotatedLibraryEntry unitEntry=model.getUnitSpectrum(filename, accessions, precursorCharge, retentionTime, params);
		
		return Optional.ofNullable(extractSpectrum(unitEntry, duration, limitToQuantifiable));
	}

	public TransitionRefinementData extractSpectrum(AnnotatedLibraryEntry unitEntry, float duration, boolean limitToQuantifiable) {
		ArrayList<Stripe> stripes=getScanSubset(unitEntry.getRetentionTime()-duration, unitEntry.getRetentionTime()+duration);
		return quantifyPeptide(scorer, unitEntry, limitToQuantifiable, stripes);
	}

	public static TransitionRefinementData quantifyPeptide(PSMScorer scorer, AnnotatedLibraryEntry unitEntry, boolean limitToQuantifiable, ArrayList<Stripe> stripes) {
		// find the center
		float bestDelta=Float.MAX_VALUE;
		PeakScores[] bestScores=null;
		ArrayList<PeakScores[]> scoreList=new ArrayList<PeakScores[]>();
		TFloatArrayList retentionTimes=new TFloatArrayList();
		TFloatArrayList totalIonCurrent=new TFloatArrayList();
		TFloatArrayList totalIdentifiedIonCurrent=new TFloatArrayList();
		for (Stripe stripe : stripes) {
			retentionTimes.add(stripe.getScanStartTime());
			float delta=Math.abs(stripe.getScanStartTime()-unitEntry.getRetentionTime());
			PeakScores[] individualPeakScores=scorer.getIndividualPeakScores(unitEntry, stripe, true);
			scoreList.add(individualPeakScores);
			if (delta<bestDelta) {
				bestDelta=delta;
				bestScores=individualPeakScores;
			}
			float sumIdentifiedIntensities=0.0f;
			for (int i=0; i<individualPeakScores.length; i++) {
				if (individualPeakScores[i]!=null) {
					sumIdentifiedIntensities+=individualPeakScores[i].getScore();
				}
			}
			totalIdentifiedIonCurrent.add(sumIdentifiedIntensities);
			totalIonCurrent.add(stripe.getTIC());
		}
		// no signal of any kind at retention time!
		if (bestScores==null||bestScores.length==0) return null;

		// get each scan (fragments by RT)
		TFloatArrayList[] traces=new TFloatArrayList[bestScores.length];
		@SuppressWarnings("unchecked")
		ArrayList<XYZPoint>[] deltaMassesByRT=new ArrayList[bestScores.length];
		for (int i=0; i<traces.length; i++) {
			traces[i]=new TFloatArrayList();
			deltaMassesByRT[i]=new ArrayList<XYZPoint>();
		}
		for (int index=0; index<scoreList.size(); index++) {
			PeakScores[] peakScores=scoreList.get(index);
			for (int i=0; i<peakScores.length; i++) {
				if (peakScores[i]!=null) {
					traces[i].add(peakScores[i].getScore());
					deltaMassesByRT[i].add(new XYZPoint(retentionTimes.get(index), peakScores[i].getDeltaMass(), peakScores[i].getScore()));
				} else {
					traces[i].add(0.0f);
				}
			}
		}

		// invert each scan into fragment chromatograms (RTs by fragment)
		ArrayList<PeakScores> bestKeptPeaks=new ArrayList<PeakScores>();
		ArrayList<float[]> chromatograms=new ArrayList<float[]>();
		ArrayList<ArrayList<XYZPoint>> chromatogramDeltaMassesByRT=new ArrayList<ArrayList<XYZPoint>>();
		ArrayList<FragmentIon> fragmentMasses=new ArrayList<FragmentIon>();
		for (int i=0; i<bestScores.length; i++) {
			if (bestScores[i]!=null&&bestScores[i].getScore()>0) {
				float[] chromatogram=traces[i].toArray();
				chromatogram=SkylineSGFilter.paddedSavitzkyGolaySmooth(chromatogram);
				chromatograms.add(chromatogram);
				chromatogramDeltaMassesByRT.add(deltaMassesByRT[i]);
				bestKeptPeaks.add(bestScores[i]);
				fragmentMasses.add(bestScores[i].getTarget());
			}
		}
		
		// identify transitions
		TransitionRefinementData data=TransitionRefiner.identifyTransitions(unitEntry.getPeptideModSeq(), unitEntry.getPrecursorCharge(), fragmentMasses.toArray(new FragmentIon[fragmentMasses.size()]), chromatograms, retentionTimes.toArray());
		float[] correlations=data.getCorrelationArray();
		float[] integrations=data.getIntegrationArray();
		Range rtRange=data.getRange();

		TDoubleArrayList mzs=new TDoubleArrayList();
		TFloatArrayList intens=new TFloatArrayList();
		TFloatArrayList deltaMasses=new TFloatArrayList(); // will ultimately be the length of the correlations array

		float correlationThreshold=limitToQuantifiable?TransitionRefiner.quantitativeCorrelationThreshold:-1f;
		for (int i=0; i<bestKeptPeaks.size(); i++) {
			// calculate delta mass for each fragment ion
			
			float totalDeltaMasses=0.0f;
			float totalIntensities=0.0f;
			for (XYZPoint point : chromatogramDeltaMassesByRT.get(i)) {
				if (rtRange.contains((float)point.getX())) {
					totalDeltaMasses+=point.getY();
					totalIntensities+=point.getZ();
				}
			}
			
			float deltaMass=0.0f;
			if (totalIntensities>0.0f) {
				deltaMass=totalDeltaMasses/totalIntensities;
			}
			deltaMasses.add(deltaMass);
			
			// generate spectrum peaks for only each "kept" fragment ion
			if (correlations[i]>=correlationThreshold) {
				// grab mz
				PeakScores bestScore=bestKeptPeaks.get(i);
				float peakScore=bestScore.getScore();
				if (peakScore>0) {
					mzs.add(bestScore.getTargetMass());
					intens.add(integrations[i]);
				}
			}
		}

		float ticSum=0.0f;
		float identifiedTicSum=0.0f;
		for (int i=0; i<stripes.size(); i++) {
			Stripe stripe=stripes.get(i);
			if (rtRange.contains(stripe.getScanStartTime())) {
				identifiedTicSum+=totalIdentifiedIonCurrent.get(i);
				ticSum+=totalIonCurrent.get(i);
			}
		}
		float identifiedTICRatio=ticSum==0.0f?0.0f:identifiedTicSum/ticSum;

		// System.out.println(peptideModSeq+"\t"+keptPeaks.size()+"\t"+count+"\t"+quantCount);

		double[] massArray=mzs.toArray();
		float[] intensityArray=intens.toArray();
		float[] deltaMassArray=deltaMasses.toArray();
		return data.addPeakData(deltaMassArray, massArray, intensityArray, retentionTimes.toArray(), identifiedTICRatio);
	}
}
