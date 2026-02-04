package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ParsimonyProteinGrouper;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.LibraryPeakLocationInferrer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorProteinGroup;
import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;
import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.NonstandardAminoAcidException;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.PeptideQuantExtractor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AbstractSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaFeaturePredictionModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaLibraryPredictionClient;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaPrecursor;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation.Prosit2019HCDModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt.Prosit2019RTModel;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class Reintegrator {
	private static final boolean WRITE_REPORTS=false; //FIXME turn off for production
	private static final boolean WRITE_LIBRARIES=true; //FIXME turn on for production
	private static final boolean INTEGRATE_PRECURSORS=false;
	
	public static void main(String[] args) throws Exception {
		File mainDir=new File("/Users/searleb/Documents/manuscripts/2025/mapms/mapms/");
		File spectronautCSV;
		File[] rawFiles;
		File newLibraryFile;
		File prositLibraryFile=new File("/Users/searleb/Documents/encyclopedia/human/uniprot_human_25apr2019.z3_nce.dlib");

		if (false) {
			spectronautCSV=new File(mainDir, "BGS_TP_SN19_125ng_1e6_6xpDIA_directDIA_Report.csv");
		 	rawFiles=new File[] {new File(mainDir, "2024_06_17_125ng_HeLa_1e6_6xpDIA_aurora_01.dia"), new File(mainDir, "2024_06_17_125ng_HeLa_1e6_6xpDIA_aurora_02.dia"), new File(mainDir, "2024_06_17_125ng_HeLa_1e6_6xpDIA_aurora_03.dia")};
			newLibraryFile=new File(mainDir, "BGS_TP_SN19_125ng_1e6_6xpDIA_directDIA_Report.elib");
			prositLibraryFile=new File(mainDir, "BGS_TP_SN19_125ng_1e6_6xpDIA_directDIA_Report.csv_predicted.dlib");
		} else if (true) {
			spectronautCSV=new File(mainDir, "BGS_TP_SN19_125ng_1e6_NormDIA_directDIA_Report.csv");
		 	rawFiles=new File[] {new File(mainDir, "2024_06_17_125ng_HeLa_1e6_NormDIA_aurora_01.dia"), new File(mainDir, "2024_06_17_125ng_HeLa_1e6_NormDIA_aurora_02.dia"), new File(mainDir, "2024_06_17_125ng_HeLa_1e6_NormDIA_aurora_03.dia")};
			newLibraryFile=new File(mainDir, "BGS_TP_SN19_125ng_1e6_NormDIA_directDIA_Report.elib");
			prositLibraryFile=new File(mainDir, "BGS_TP_SN19_125ng_1e6_NormDIA_directDIA_Report.csv_predicted.dlib");
		} else {
			spectronautCSV=new File(mainDir, "BGS_TP_SN19_125ng_3e6_NormDIA_directDIA_Report.csv");
		 	rawFiles=new File[] {new File(mainDir, "2024_06_17_125ng_HeLa_3e6_NormDIA_aurora_01.dia"), new File(mainDir, "2024_06_17_125ng_HeLa_3e6_NormDIA_aurora_02.dia"), new File(mainDir, "2024_06_17_125ng_HeLa_3e6_NormDIA_aurora_03.dia")};
			newLibraryFile=new File(mainDir, "BGS_TP_SN19_125ng_3e6_NormDIA_directDIA_Report.elib");
			prositLibraryFile=new File(mainDir, "BGS_TP_SN19_125ng_3e6_NormDIA_directDIA_Report.csv_predicted.dlib");
		}
		
		HashMap<String, String> defaultParameters = SearchParameterParser.getDefaultParameters();
		defaultParameters.put(SearchParameters.SMOOTH_INTEGRATIONS, Boolean.FALSE.toString());
		defaultParameters.put(SearchParameters.SUBTRACT_BACKGROUND, Boolean.FALSE.toString());
		defaultParameters.put("-ptol", 5+"");
		final SearchParameters params=SearchParameterParser.parseParameters(defaultParameters);
		
		reintegrateFromSpectronaut(spectronautCSV, rawFiles, newLibraryFile, Optional.of(prositLibraryFile), params);
	}
	
	public static class Dataset {
		final HashMap<String, HashSet<String>> targetAccessionsByPeptide=new HashMap<>();
		final HashMap<String, File> filesByBasename=new HashMap<String, File>();
		final HashMap<String, ArrayList<PSMData>> matchesMap=new HashMap<String, ArrayList<PSMData>>();
		final ArrayList<PercolatorPeptide> scoredPeptides=new ArrayList<PercolatorPeptide>();
		final ArrayList<PercolatorPeptide> decoyPeptides=new ArrayList<PercolatorPeptide>();
	}

	public static void reintegrateFromSpectronaut(final File spectronautCSV, final File[] rawFiles, File newLibraryFile, Optional<File> prositLibraryFile, SearchParameters params)
			throws IOException, SQLException, DataFormatException, InterruptedException {
		Pair<Dataset, LibraryPeakLocationInferrer> result = readDataset(spectronautCSV, rawFiles, newLibraryFile, prositLibraryFile, params);

		LibraryFile newLibrary=null;
		if (WRITE_LIBRARIES) {
			newLibrary=new LibraryFile();
			newLibrary.openFile();
			Logger.logLine("Writing library file "+newLibrary.getName());
			newLibrary.dropIndices();
		}
		
		TFloatArrayList[] intensityValuesByCorrelation=new TFloatArrayList[101];
		for (int i = 0; i < intensityValuesByCorrelation.length; i++) {
			intensityValuesByCorrelation[i]=new TFloatArrayList();
		}

		int[] pointsAcrossThePeakHistogram=new int[21];
		HashMap<String, TFloatArrayList> quantitativeValuesByPeptide=new HashMap<String, TFloatArrayList>();
		
		ProgressIndicator progress=new EmptyProgressIndicator();
		for (Entry<String, ArrayList<PSMData>> entryMapping : result.x.matchesMap.entrySet()) {
			ArrayList<PSMData> matches=entryMapping.getValue();
			File rawFile=result.x.filesByBasename.get(entryMapping.getKey());
			
			StripeFileInterface stripeFile=StripeFileGenerator.getFile(rawFile, params, true);
			
			PeptideQuantExtractor extractor=new PeptideQuantExtractor(progress, stripeFile, params);
			ArrayList<IntegratedLibraryEntry> entries=extractor.extractPeptides(matches, Optional.of(result.y), false);
			
			for (IntegratedLibraryEntry entry : entries) {
				String peptideKey=getPeptideKey(entry.getPeptideModSeq(), entry.getPrecursorCharge());
				Optional<Float> correlationWithFragments = entry.getRefinementData().getCorrelationWithFragments();
				if (correlationWithFragments.isPresent()) {
					int index=Math.round(100f*correlationWithFragments.get());
					if (index<0) index=0;
					if (index>100) index=100;
					Float totalIntensity = entry.getRefinementData().getTotalPrecursorIntensity().get();
					intensityValuesByCorrelation[index].add(totalIntensity);
					
					Optional<ArrayList<XYPoint>[]> chromatograms = entry.getRefinementData().getPrecursorChromatograms();
					ArrayList<XYPoint> monoisotopic=chromatograms.get()[1];
					
					float[] intensities=XYTrace.toFloatArrays(monoisotopic).y;
					float max=0;
					int maxIndex=0;
					for (int i = 0; i < intensities.length; i++) {
						if (intensities[i]>max) {
							max=intensities[i];
							maxIndex=i;
						}
					}
					
					int pointsAcrossThePeak=1;
					
					for (int i = maxIndex+1; i < intensities.length; i++) {
						if (intensities[i]>0) {
							pointsAcrossThePeak++;
						} else {
							break;
						}
					}
					for (int i = maxIndex-1; i >= 0; i--) {
						if (intensities[i]>0) {
							pointsAcrossThePeak++;
						} else {
							break;
						}
					}
					
					if (pointsAcrossThePeak>=pointsAcrossThePeakHistogram.length) {
						pointsAcrossThePeak=pointsAcrossThePeakHistogram.length-1;
					}
					pointsAcrossThePeakHistogram[pointsAcrossThePeak]++;
					
					TFloatArrayList quantValues=quantitativeValuesByPeptide.get(peptideKey);
					if (quantValues==null) {
						quantValues=new TFloatArrayList();
						quantitativeValuesByPeptide.put(peptideKey, quantValues);
					}
					quantValues.add(totalIntensity);
				}
			}

			if (WRITE_LIBRARIES) {
				newLibrary.addTIC(stripeFile);
				
				newLibrary.addIntegratedEntries(!INTEGRATE_PRECURSORS, entries, Optional.empty(), Optional.empty(), params);		
			}
		}

		if (WRITE_REPORTS) {
			System.out.println("correlation\tcount\t5p\t25p\t50p\t75p\t95p");
			for (int i = 0; i < intensityValuesByCorrelation.length; i++) {
				float[] intensities=intensityValuesByCorrelation[i].toArray();
	
				System.out.println(i + "\t" + intensities.length + "\t" + QuickMedian.select(intensities, 0.05f) + "\t"
						+ QuickMedian.select(intensities, 0.25f) + "\t" + QuickMedian.select(intensities, 0.5f) + "\t"
						+ QuickMedian.select(intensities, 0.75f) + "\t" + QuickMedian.select(intensities, 0.95f));
			}
			System.out.println();
			
			System.out.println("points\tcount");
			for (int i = 0; i < pointsAcrossThePeakHistogram.length; i++) {
				System.out.println(i+"\t"+pointsAcrossThePeakHistogram[i]);
			}
			System.out.println();
			
			float[] cvs=new float[300];
			for (TFloatArrayList list : quantitativeValuesByPeptide.values()) {
				float[] data=list.toArray();
				float stdev = General.stdev(data);
				if (stdev>0) {
					float mean = General.mean(data);
					float cv=stdev/mean;
					int index=Math.round(100f*cv);
					if (index>=cvs.length) index=cvs.length-1;
					cvs[index]++;
				}
			}
			
			System.out.println("cv\tcount");
			for (int i = 0; i < cvs.length; i++) {
				System.out.println((1*i)+"\t"+cvs[i]);
			}
			System.out.println();
		}
		
		PrintWriter writer=new PrintWriter(new File(newLibraryFile.getParentFile(), spectronautCSV.getName()+"_averages.txt"));
		for (Entry<String, TFloatArrayList> entry : quantitativeValuesByPeptide.entrySet()) {
			float average = General.sum(entry.getValue().toArray())/rawFiles.length;
			float loggedValue=average<=0?0:Log.log10(average); 
			writer.println(entry.getKey()+"\t"+loggedValue);
		}
		writer.flush();
		writer.close();

		if (WRITE_LIBRARIES) {
			newLibrary.addProteinsFromEntries(result.x.targetAccessionsByPeptide, new HashMap<String, HashSet<String>>());
			newLibrary.addTargetDecoyPeptides(result.x.scoredPeptides, result.x.decoyPeptides);
			Pair<ArrayList<PercolatorProteinGroup>, ArrayList<PercolatorProteinGroup>> targetDecoyProteins=ParsimonyProteinGrouper.groupProteins(result.x.scoredPeptides, result.x.decoyPeptides, params.getPercolatorProteinThreshold(), params.getAAConstants());
	
			Logger.logLine("Writing global target/decoy proteins: "+targetDecoyProteins.x.size()+"/"+targetDecoyProteins.y.size());
			newLibrary.addTargetDecoyProteins(SearchToBLIB.GLOBAL_NAME, targetDecoyProteins.x, targetDecoyProteins.y);
	
			newLibrary.createIndices();
			newLibrary.saveAsFile(newLibraryFile);
		}
	}
	
	private static class SpectronautJobData extends AbstractSearchJobData {
		private SpectronautJobData(File diaFile, SearchParameters params) {
			super(diaFile, null, params, ProgramType.getGlobalVersion().toString());
		}

		@Override
		public String getSearchType() {
			return "Spectronaut";
		}

		@Override
		public String getPrimaryScoreName() {
			return "EG.Qvalue";
		}

		@Override
		public SearchJobData updateQuantFile(File f) {
			return new SpectronautJobData(f, getParameters());
		}
		
		
	}

	public static Pair<Dataset, LibraryPeakLocationInferrer> readDataset(final File spectronautCSV,
			final File[] rawFiles, File newLibraryFile, Optional<File> prositLibraryFile, SearchParameters params)
			throws IOException, SQLException, DataFormatException {
		Dataset dataset=new Dataset();
		
		HashMap<String, SpectronautJobData> jobDataByFile=new HashMap<String, Reintegrator.SpectronautJobData>();
		HashMap<SearchJobData, TObjectFloatHashMap<String>> rtByPeptideModSeq=new HashMap<SearchJobData, TObjectFloatHashMap<String>>();
		
		for (int i = 0; i < rawFiles.length; i++) {
			String basename=rawFiles[i].getName();
			basename=basename.substring(0, basename.lastIndexOf('.'));
			dataset.filesByBasename.put(basename, rawFiles[i]);
			dataset.matchesMap.put(basename, new ArrayList<PSMData>());
			
			SpectronautJobData jobData=new SpectronautJobData(rawFiles[i], params);
			jobDataByFile.put(basename, jobData);
			rtByPeptideModSeq.put(jobData, new TObjectFloatHashMap<String>());
		}
		
		
		Logger.logLine("Reading existing predictions from library");
		final HashMap<String, LibraryEntry> previousPredictions=new HashMap<String, LibraryEntry>();
		if (prositLibraryFile.isPresent()) {
			LibraryFile prositLibrary=new LibraryFile();
			prositLibrary.openFile(prositLibraryFile.get());
			ArrayList<LibraryEntry> entries=prositLibrary.getAllEntries(false, params.getAAConstants());
			for (LibraryEntry entry : entries) {
				previousPredictions.put(getPeptideKey(entry.getPeptideModSeq(), entry.getPrecursorCharge()), entry);
			}
		}
		Logger.logLine("Found "+previousPredictions.size()+" total previous predictions");
		
		final ArrayList<PSMData> requiresPrediction=new ArrayList<PSMData>();
		final HashMap<String, LibraryEntry> selectedPredictions=new HashMap<String, LibraryEntry>();
		
		TableParser.parseCSV(spectronautCSV, new TableParserMuscle() {
			
			@Override
			public void processRow(Map<String, String> row) {
				String filename=row.get("R.FileName");
				
				ArrayList<PSMData> matches=dataset.matchesMap.get(filename);
				if (matches==null) {
					return;
				}
				File rawFile=dataset.filesByBasename.get(filename);
				
				String accessions=row.get("PG.ProteinAccessions");
				HashSet<String> accessionSet=new HashSet<String>(Arrays.asList(accessions.split(";")));
				String intPIMID=row.get("FG.IntMID");
				intPIMID=intPIMID.substring(1, intPIMID.length() - 1);
				String peptideModSeq=PeptideUtils.getCorrectedMasses(intPIMID);
				
				String mzString=row.get("FG.PrecMz");
				double precursorMZ=Double.parseDouble(mzString);
				
				String chargeString=row.get("FG.Charge");
				byte precursorCharge=Byte.parseByte(chargeString);
				
				String rtString=row.get("EG.ApexRT");
				float retentionTime=Float.parseFloat(rtString)*60f; // in sec
				
				String scoreString=row.get("FG.CScore");
				float score=Float.parseFloat(scoreString);
				
				String PEPString=row.get("EG.PEP");
				float PEP=Float.parseFloat(PEPString);
				
				String qvalueString=row.get("EG.Qvalue");
				float qValue=Float.parseFloat(qvalueString);
				
				int spectrumIndex=0;
				float duration = params.getExpectedPeakWidth();

				PSMData match = new PSMData(accessionSet, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq,
						retentionTime, score, score, duration, false, params.getAAConstants());
				
				LibraryEntry prediction=previousPredictions.get(getPeptideKey(peptideModSeq, precursorCharge));
				if (prediction!=null) {
					selectedPredictions.put(getPeptideKey(peptideModSeq, precursorCharge), prediction);
				} else {
					requiresPrediction.add(match);
				}
				
				matches.add(match);
				
				dataset.targetAccessionsByPeptide.put(match.getPeptideSeq(), match.getAccessions());
				
				String psmID=PercolatorPeptide.getPSMID(rawFile.getName(), retentionTime, Optional.empty(), false, peptideModSeq, precursorCharge);
				PercolatorPeptide pep=new PercolatorPeptide(psmID, PSMData.accessionsToString(accessionSet), 
						qValue, PEP, params.getAAConstants());
				dataset.scoredPeptides.add(pep);
				
				SpectronautJobData jobData=jobDataByFile.get(filename);
				TObjectFloatHashMap<String> rtMap=rtByPeptideModSeq.get(jobData);
				rtMap.put(peptideModSeq, retentionTime);
			}
			
			@Override
			public void cleanup() {
			}
		});
		
		for (Entry<String, ArrayList<PSMData>> entryMapping : dataset.matchesMap.entrySet()) {
			Logger.logLine("Found "+entryMapping.getValue().size()+" matches from "+entryMapping.getKey());
		}
		
		Logger.logLine("Matches correspond to "+dataset.targetAccessionsByPeptide.size()+" unique peptides");
		
		Logger.logLine("Found "+selectedPredictions.size()+" peptides matching existing predictions. Need to add "+requiresPrediction.size()+" new predictions");

		ArrayList<KoinaFeaturePredictionModel> models=new ArrayList<KoinaFeaturePredictionModel>();
		models.add(new Prosit2019HCDModel());
		models.add(new Prosit2019RTModel());
		KoinaLibraryPredictionClient client=new KoinaLibraryPredictionClient(models);
		
		ArrayList<KoinaPrecursor> pred=new ArrayList<KoinaPrecursor>();
		for (PSMData data : requiresPrediction) {
			try {
				AminoAcidEncoding[] encoding = AminoAcidEncoding.getAAs(data.getPeptideModSeq(), params.getAAConstants());
				boolean passes=client.checkPeptide(encoding, data.getPrecursorCharge());
			
				KoinaPrecursor precursor = new KoinaPrecursor(encoding, 33f, data.getPrecursorCharge());
				if (passes) {
					pred.add(precursor);
				}
			} catch (NonstandardAminoAcidException e) {
				Logger.errorLine("FAILED: "+data.getPeptideModSeq());
			}
		}
		client.generatePredictions(KoinaLibraryPredictionClient.HTTPS_KOINA_WILHELMLAB_ORG_443, pred, new EmptyProgressIndicator(true));
		
		
		
		for (KoinaPrecursor precursor : pred) {
			AnnotatedLibraryEntry entry=precursor.toEntry(params.getAAConstants(), params);
			selectedPredictions.put(getPeptideKey(entry.getPeptideModSeq(), entry.getPrecursorCharge()), entry);
		}
		
		// FIXME THIS INFERRER IS NOT SET UP!!!!
		LibraryPeakLocationInferrer inferrer=new LibraryPeakLocationInferrer(selectedPredictions.values(), rtByPeptideModSeq, params);
		Logger.logLine("Finished predictions! Able to capture "+selectedPredictions.size()+" total predictions");

		if (false) {
			// NOTE: this section previously used for setting up libraries for faster analysis
			LibraryFile predicted=new LibraryFile();
			predicted.openFile();
			Logger.logLine("Writing library file "+predicted.getName());
			predicted.dropIndices();
			ArrayList<LibraryEntry> returnList = new ArrayList<LibraryEntry>(selectedPredictions.values());
			predicted.addEntries(returnList, false);
			predicted.addProteinsFromEntries(returnList);
			predicted.createIndices();
			predicted.saveAsFile(new File(newLibraryFile.getParentFile(), spectronautCSV.getName()+"_predicted.dlib"));
		}
		Pair<Dataset, LibraryPeakLocationInferrer> result=new Pair<Reintegrator.Dataset, LibraryPeakLocationInferrer>(dataset, inferrer);
		return result;
	}
	
	private static String getPeptideKey(String peptideModSeq, byte precursorCharge) {
		return peptideModSeq+"_"+precursorCharge;
	}
}
