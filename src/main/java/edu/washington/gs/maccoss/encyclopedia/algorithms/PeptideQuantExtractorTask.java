package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Optional;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoPermuter;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SkylineSGFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ThreadableTask;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.set.hash.TFloatHashSet;

public class PeptideQuantExtractorTask extends ThreadableTask<Nothing> {
	private final Optional<LibraryInterface> library;
	private final ArrayList<Stripe> stripes;
	private final boolean limitToQuantifiable;

	private final PSMScorer scorer;
	private final SearchParameters params;

	private final PSMData psmdata;
	private final ConcurrentLinkedQueue<IntegratedLibraryEntry> savedEntries; // CAN BE NULL

	public PeptideQuantExtractorTask(PSMData psmdata, Optional<LibraryInterface> library, ArrayList<Stripe> stripes, SearchParameters parameters, boolean limitToQuantifiable) {
		this.psmdata=psmdata;
		this.library=library;
		this.stripes=stripes;

		scorer=new DotProduct(parameters.getFragmentTolerance());
		params=parameters;
		this.savedEntries=null;
		
		this.limitToQuantifiable=limitToQuantifiable; //library.isPresent();
	}

	public PeptideQuantExtractorTask(PSMData psmdata, Optional<LibraryInterface> library, ArrayList<Stripe> stripes, SearchParameters parameters, ConcurrentLinkedQueue<IntegratedLibraryEntry> savedEntries, boolean limitToQuantifiable) {
		this.psmdata=psmdata;
		this.library=library;
		this.stripes=stripes;

		scorer=new DotProduct(parameters.getFragmentTolerance());
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
		Optional<TransitionRefinementData> spectrum=extractSpectrum(psmdata.getAccessions(), library, psmdata.getPrecursorCharge(), psmdata.getPeptideModSeq(), psmdata.getRetentionTime(), psmdata.getDuration(), limitToQuantifiable);
		if (params.isRunPhosphoLocalization()) {
			ArrayList<String> permutations=PhosphoPermuter.getPermutations(psmdata.getPeptideModSeq(), params.getAAConstants());
			if (permutations.size()==1) {
				System.out.println("single\t"+psmdata.getPeptideModSeq());
			} else {
				boolean multiple=extractPhosphoForms(psmdata.getPrecursorMZ(), psmdata.getPrecursorCharge(), permutations, psmdata.getRetentionTime());
				System.out.println("multiple\t"+psmdata.getPeptideModSeq()+"\t"+multiple);
			}
		}
		if (spectrum.isPresent()) {
			// FIXME need to not add duplicates!!!! for now just run SQL:
			// delete from entries where RowId not in (SELECT MIN(RowId) FROM entries GROUP BY PeptideModSeq, PrecursorCharge)
			TransitionRefinementData data=spectrum.get();
			IntegratedLibraryEntry entry=new IntegratedLibraryEntry(psmdata.getAccessions(), psmdata.getSpectrumIndex(), psmdata.getPrecursorMZ(), psmdata.getPrecursorCharge(), psmdata.getPeptideModSeq(), 1, psmdata.getRetentionTime(), psmdata.getScore(), data.getMassArray().get(), data.getIntensityArray().get(), data.getRange());
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



	private boolean extractPhosphoForms(double precursorMZ, byte precursorCharge, ArrayList<String> peptideModSeqs, float retentionTime) {
		float duration=6*60f; // search for 6 minutes
		
		EncyclopediaOneScorer encyclopediaScorer=new EncyclopediaOneScorer(params, null); // not using aux scoring

		ArrayList<Stripe> stripes=getScanSubset(retentionTime-duration, retentionTime+duration);

		TFloatArrayList allBestTimes=new TFloatArrayList();

		TFloatHashSet list=new TFloatHashSet();
		for (String peptideModSeq : peptideModSeqs) {
			FragmentationModel model=new FragmentationModel(peptideModSeq, params.getAAConstants());
			LibraryEntry unitEntry=model.getUnitSpectrum(new HashSet<String>(), precursorCharge, retentionTime, params);

			TFloatArrayList bestTimes=new TFloatArrayList();
			float bestScore=0.0f;
			TFloatFloatHashMap rtScoreMap=new TFloatFloatHashMap();
			for (Stripe spectrum : stripes) {
				// float score=encyclopediaScorer.score(unitEntry, spectrum);
				PeakScores[] individualPeakScores=encyclopediaScorer.getIndividualPeakScores(unitEntry, spectrum, false);
				float score=0.0f;
				for (int i=0; i<individualPeakScores.length; i++) {
					if (individualPeakScores[i]!=null) {
						score++;
					}
				}
				rtScoreMap.put(spectrum.getScanStartTime(), score);
				if (score>9.5f) { // at least 5 transitions
					if (bestScore<score) {
						bestScore=score;
						bestTimes.clear();
						bestTimes.add(spectrum.getScanStartTime());
					} else if (bestScore==score) {
						bestTimes.add(spectrum.getScanStartTime());
					}
				}
			}

			allBestTimes.addAll(bestTimes);

			EValueCalculator calculator=new EValueCalculator(rtScoreMap);
			// System.out.println(peptideModSeq+"\t"+calculator.getMaxRT()+"\t"+calculator.getNegLog10EValue()+"\t"+calculator.getMaxRawScore()); //FIXME
			list.add(calculator.getMaxRT());
		}

		if (allBestTimes.size()>1) {
			float range=allBestTimes.max()-allBestTimes.min();
			if (range>=60) {
				// System.out.println("multiple\t"+peptideModSeqs.get(0));
				return true;
			}
		} else if (allBestTimes.size()==1) {
			// System.out.println("single\t"+peptideModSeqs.get(0));
		}
		return false;
	}

	private Optional<TransitionRefinementData> extractSpectrum(HashSet<String> accessions, Optional<LibraryInterface> library, byte precursorCharge, String peptideModSeq, float retentionTime, float duration, boolean limitToQuantifiable) {
		LibraryEntry unitEntry=null;
		if (library.isPresent()) {
			try {
				ArrayList<LibraryEntry> entries=library.get().getEntries(peptideModSeq, precursorCharge, false);
				if (entries.size()>0) {
					unitEntry=entries.get(0).toUnitSpectrum();
				} else {
					 // if library is ok but spectrum is not in library, just return null (don't quantify)
					return Optional.empty();
				}
				if (unitEntry.getIonCount()<1) { // FIXME HACKS GALORE
					 // if unit spectrum has fewer than 1 peaks, just return null (don't quantify)
					return Optional.empty();
				}
			} catch (IOException ioe) {
				Logger.errorLine("Error processing "+library.get().getName());
				throw new EncyclopediaException("Error parsing Stripe file", ioe);
			} catch (SQLException sqle) {
				Logger.errorLine("Error processing "+library.get().getName());
				throw new EncyclopediaException("Error parsing Stripe file", sqle);
			} catch (DataFormatException dfe) {
				Logger.errorLine("Error processing "+library.get().getName());
				throw new EncyclopediaException("Error parsing Stripe file", dfe);
			}
		}
		
		if (unitEntry==null) {
			FragmentationModel model=new FragmentationModel(peptideModSeq, params.getAAConstants());
			unitEntry=model.getUnitSpectrum(accessions, precursorCharge, retentionTime, params);
		}
		
		return Optional.ofNullable(extractSpectrum(unitEntry, duration, limitToQuantifiable));
	}

	public TransitionRefinementData extractSpectrum(LibraryEntry unitEntry, float duration, boolean limitToQuantifiable) {
		ArrayList<Stripe> stripes=getScanSubset(unitEntry.getRetentionTime()-duration, unitEntry.getRetentionTime()+duration);
		return quantifyPeptide(scorer, unitEntry, limitToQuantifiable, stripes);
	}

	public static TransitionRefinementData quantifyPeptide(PSMScorer scorer, LibraryEntry unitEntry, boolean limitToQuantifiable, ArrayList<Stripe> stripes) {
		float bestDelta=Float.MAX_VALUE;
		PeakScores[] bestScores=null;
		ArrayList<PeakScores[]> scoreList=new ArrayList<PeakScores[]>();
		TFloatArrayList retentionTimes=new TFloatArrayList();
		for (Stripe stripe : stripes) {
			retentionTimes.add(stripe.getScanStartTime());
			float delta=Math.abs(stripe.getScanStartTime()-unitEntry.getRetentionTime());
			PeakScores[] individualPeakScores=scorer.getIndividualPeakScores(unitEntry, stripe, true);
			scoreList.add(individualPeakScores);
			if (delta<bestDelta) {
				bestDelta=delta;
				bestScores=individualPeakScores;
			}
		}
		// no signal of any kind at retention time!
		if (bestScores==null) return null;

		TFloatArrayList[] traces=new TFloatArrayList[bestScores.length];
		for (int i=0; i<traces.length; i++) {
			traces[i]=new TFloatArrayList();
		}
		for (PeakScores[] peakScores : scoreList) {
			for (int i=0; i<peakScores.length; i++) {
				if (peakScores[i]!=null) {
					traces[i].add(peakScores[i].getScore());
				} else {
					traces[i].add(0.0f);
				}
			}
		}

		ArrayList<PeakScores> keptPeaks=new ArrayList<PeakScores>();
		ArrayList<float[]> chromatograms=new ArrayList<float[]>();
		for (int i=0; i<bestScores.length; i++) {
			if (bestScores[i]!=null&&bestScores[i].getScore()>0) {
				float[] chromatogram=traces[i].toArray();
				chromatogram=SkylineSGFilter.paddedSavitzkyGolaySmooth(chromatogram);
				chromatograms.add(chromatogram);
				keptPeaks.add(bestScores[i]);
			}
		}
		
		TransitionRefinementData data=TransitionRefiner.identifyTransitions(unitEntry.getPeptideModSeq(), chromatograms, retentionTimes.toArray());
		float[] correlations=data.getCorrelationArray();
		float[] integrations=data.getIntegrationArray();

		TDoubleArrayList mzs=new TDoubleArrayList();
		TFloatArrayList intens=new TFloatArrayList();
		//int count=0;
		//int quantCount=0;
		float correlationThreshold=limitToQuantifiable?TransitionRefiner.quantitativeCorrelationThreshold:0.0f;//TransitionRefiner.identificationCorrelationThreshold;
		for (int i=0; i<keptPeaks.size(); i++) {
			PeakScores scores=keptPeaks.get(i);
			if (correlations[i]>=correlationThreshold) {
				//quantCount++;
				float peakScore=scores.getScore();
				if (peakScore>0) {
					mzs.add(scores.getTargetMass());
					intens.add(integrations[i]);
					//count++;
				}
			}
		}
		
		if (mzs.size()==0) return null;

		// System.out.println(peptideModSeq+"\t"+keptPeaks.size()+"\t"+count+"\t"+quantCount);

		double[] massArray=mzs.toArray();
		float[] intensityArray=intens.toArray();
		return data.addPeakData(massArray, intensityArray, retentionTimes.toArray());
	}
}
