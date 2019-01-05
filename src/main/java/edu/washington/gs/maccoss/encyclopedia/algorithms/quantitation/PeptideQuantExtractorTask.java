package edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Optional;
import java.util.concurrent.ConcurrentLinkedQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.DotProduct;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMPeakScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeakLocationInferrerInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.AmbiguousPeptideModSeq;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.*;
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
	private final Optional<PeakLocationInferrerInterface> inferrer;
	private final String filename;
	private final ArrayList<Stripe> stripes;
	private final boolean limitToQuantifiable;

	private final PSMPeakScorer scorer;
	private final SearchParameters params;

	private final PSMData psmdata;
	private final ConcurrentLinkedQueue<IntegratedLibraryEntry> savedEntries; // CAN BE NULL

	public PeptideQuantExtractorTask(String filename, PSMData psmdata, final Optional<PeakLocationInferrerInterface> inferrer, Optional<PhosphoLocalizer> localizer, ArrayList<Stripe> stripes, SearchParameters parameters, boolean limitToQuantifiable) {
		this.filename=filename;
		this.psmdata=psmdata;
		this.inferrer=inferrer;
		this.localizer=localizer;
		this.stripes=stripes;

		scorer=new DotProduct(parameters);
		params=parameters;
		this.savedEntries=null;
		
		this.limitToQuantifiable=limitToQuantifiable; //library.isPresent();
	}

	public PeptideQuantExtractorTask(String filename, PSMData psmdata, final Optional<PeakLocationInferrerInterface> inferrer, Optional<PhosphoLocalizer> localizer, ArrayList<Stripe> stripes, SearchParameters parameters, ConcurrentLinkedQueue<IntegratedLibraryEntry> savedEntries, boolean limitToQuantifiable) {
		this.filename=filename;
		this.psmdata=psmdata;
		this.inferrer=inferrer;
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
		Optional<TransitionRefinementData> spectrum=extractSpectrum(psmdata.getAccessions(), psmdata.getPrecursorCharge(), psmdata.getPeptideModSeq(), psmdata.getRetentionTime(), psmdata.getDuration(), limitToQuantifiable, inferrer, params.isQuantifySameFragmentsAcrossSamples(), psmdata.wasInferred());
		Optional<HashMap<String, TransitionRefinementData>> phosphoData=Optional.empty();
		if (canRunLocalization()) {
			Optional<PhosphoLocalizationData> localizationData=runLocalization(false);
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
					int numberOfMods=PeptideUtils.getNumberOfMods(psmdata.getPeptideModSeq(), params.getLocalizingModification().get().getNominalMass());
					if (numberOfMods>0) {
						FragmentIon[] localizingIons=new FragmentIon[0];
						float localizationScore=1000.0f;
						float numIdentificationPeaks=0.0f;
						data.setModificationLocalizationData(Optional.of(new ModificationLocalizationData(AmbiguousPeptideModSeq.getUnambigous(psmdata.getPeptideModSeq(), params.getLocalizingModification().get(), params.getAAConstants(), ""), data.getApexRT(), localizationScore, numIdentificationPeaks, numberOfMods, true, true, false, localizingIons, 0.0f, 0.0f)));
					}
				} else if (phosphoData.get().size()==0) {
					// no confident localizations, so annotate this directly on the data object
					int numberOfMods=PeptideUtils.getNumberOfMods(psmdata.getPeptideModSeq(), params.getLocalizingModification().get().getNominalMass());
					if (numberOfMods>0) {
						FragmentIon[] localizingIons=new FragmentIon[0];
						float localizationScore=0.0f;
						float numIdentificationPeaks=0.0f;
						data.setModificationLocalizationData(Optional.of(new ModificationLocalizationData(AmbiguousPeptideModSeq.getFullyAmbiguous(psmdata.getPeptideModSeq(), params.getLocalizingModification().get(), params.getAAConstants(), ""), data.getApexRT(), localizationScore, numIdentificationPeaks, numberOfMods, false, false, true, localizingIons, 0.0f, 0.0f)));
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
		return params.getLocalizingModification().isPresent()&&localizer.isPresent();
	}

	public Optional<PhosphoLocalizationData> runLocalization(boolean tryAllPermutations) {
		return localizer.get().runDIAPhosphoLocalization(psmdata, stripes, tryAllPermutations, true);
	}

	private Optional<TransitionRefinementData> extractSpectrum(HashSet<String> accessions, byte precursorCharge, String peptideModSeq, float retentionTime, float duration, boolean limitToQuantifiable, final Optional<PeakLocationInferrerInterface> inferrer, boolean integrateEverything, boolean wasInferred) {
		FragmentationModel model=PeptideUtils.getPeptideModel(peptideModSeq, params.getAAConstants());
		
		// if inferrer is present then we need to integrate everything in the target mass list
		double[] masses=null; // getUnitSpectrum is null tolerant!
		if (inferrer.isPresent()) {
			masses=inferrer.get().getTopNBestIons(peptideModSeq, precursorCharge);
		}
		AnnotatedLibraryEntry unitEntry=model.getUnitSpectrum(filename, accessions, precursorCharge, retentionTime, params, masses, 0.0, false, true);
		
		return Optional.ofNullable(extractSpectrum(unitEntry, duration, limitToQuantifiable, integrateEverything, wasInferred));
	}

	public TransitionRefinementData extractSpectrum(AnnotatedLibraryEntry unitEntry, float duration, boolean limitToQuantifiable, boolean integrateEverything, boolean wasInferred) {
		// widened to 3x the expected size (+/-1.5) of the peak to make sure we don't miss something by catching just the corner 
		ArrayList<Stripe> stripes=getScanSubset(unitEntry.getRetentionTime()-duration*1.5f, unitEntry.getRetentionTime()+duration*1.5f);
		return quantifyPeptide(scorer, unitEntry, limitToQuantifiable, stripes, integrateEverything, wasInferred, params.getAAConstants());
	}

	public static TransitionRefinementData quantifyPeptide(PSMPeakScorer scorer, AnnotatedLibraryEntry unitEntry, boolean limitToQuantifiable, ArrayList<Stripe> stripes, boolean integrateEverything, boolean wasInferred, AminoAcidConstants aaConstants) {
		// find the center
		float bestDelta=Float.MAX_VALUE;
		PeakScores[] bestScores=null;
		PeakScores[] bestIndividualFragmentScores=null;
		ArrayList<PeakScores[]> scoreList=new ArrayList<PeakScores[]>();
		TFloatArrayList retentionTimes=new TFloatArrayList();
		TFloatArrayList totalIonCurrent=new TFloatArrayList();
		TFloatArrayList totalIdentifiedIonCurrent=new TFloatArrayList();
		for (Stripe stripe : stripes) {
			retentionTimes.add(stripe.getScanStartTime());
			float delta=Math.abs(stripe.getScanStartTime()-unitEntry.getRetentionTime());
			PeakScores[] eachPeakScores=scorer.getIndividualPeakScores(unitEntry, stripe, false);
			scoreList.add(eachPeakScores);
			if (delta<bestDelta) {
				bestDelta=delta;
				bestScores=eachPeakScores;
			}
			if (bestIndividualFragmentScores==null) {
				bestIndividualFragmentScores=new PeakScores[eachPeakScores.length];
			}
			for (int i=0; i<eachPeakScores.length; i++) {
				if (bestIndividualFragmentScores[i]==null) {
					if (eachPeakScores[i]!=null) {
						bestIndividualFragmentScores[i]=eachPeakScores[i];
					}
				} else if (eachPeakScores[i]!=null&&eachPeakScores[i].getScore()>bestIndividualFragmentScores[i].getScore()) {
					bestIndividualFragmentScores[i]=eachPeakScores[i];
				}
			}
			float sumIdentifiedIntensities=0.0f;
			for (int i=0; i<eachPeakScores.length; i++) {
				if (eachPeakScores[i]!=null) {
					sumIdentifiedIntensities+=eachPeakScores[i].getScore();
				}
			}
			totalIdentifiedIonCurrent.add(sumIdentifiedIntensities);
			totalIonCurrent.add(stripe.getTIC());
			
		}
		if (integrateEverything) {
			// bestIndividualScores are all at different retention times, but all peaks are present if there is any signal for that transition
			bestScores=bestIndividualFragmentScores; 
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
		TransitionRefinementData data=TransitionRefiner.identifyTransitions(unitEntry.getPeptideModSeq(), unitEntry.getPrecursorCharge(), unitEntry.getScanStartTime(), fragmentMasses.toArray(new FragmentIon[fragmentMasses.size()]), chromatograms, retentionTimes.toArray(), wasInferred, aaConstants);
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
