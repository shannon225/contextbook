package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ParsimonyProteinGrouper;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.LibraryPeakLocationInferrer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorProteinGroup;
import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;
import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.NonstandardAminoAcidException;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.PeptideQuantExtractor;
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
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class Reintegrator {
	public static void main(String[] args) throws Exception {
		final File spectronautCSV=new File("/Users/searleb/Documents/manuscripts/2025/mapms/mapms/BGS_TP_SN19_125ng_1e6_6xpDIA_directDIA_Report.csv");
		final File rawFile=new File("/Users/searleb/Documents/manuscripts/2025/mapms/mapms/2024_06_17_125ng_HeLa_1e6_6xpDIA_aurora_01.dia");
		final File prositLibraryFile=new File("/Users/searleb/Documents/encyclopedia/human/uniprot_human_25apr2019.z3_nce.dlib");
		
		File newLibraryFile=new File("/Users/searleb/Documents/manuscripts/2025/mapms/mapms/BGS_TP_SN19_125ng_1e6_6xpDIA_directDIA_Report.elib");
		
		reintegrateFromSpectronaut(spectronautCSV, rawFile, newLibraryFile, Optional.of(prositLibraryFile));
	}

	public static void reintegrateFromSpectronaut(final File spectronautCSV, final File rawFile, File newLibraryFile, Optional<File> prositLibraryFile)
			throws IOException, SQLException, DataFormatException, InterruptedException {
		boolean integratePrecursors=false;
		ProgressIndicator progress=new EmptyProgressIndicator();
		final SearchParameters params=SearchParameterParser.getDefaultParametersObject();

		final HashMap<String, HashSet<String>> targetAccessionsByPeptide=new HashMap<>();
		final ArrayList<PSMData> matches=new ArrayList<PSMData>();
		final ArrayList<PercolatorPeptide> scoredPeptides=new ArrayList<PercolatorPeptide>();
		final ArrayList<PercolatorPeptide> decoyPeptides=new ArrayList<PercolatorPeptide>();
		
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
				if (!rawFile.getName().startsWith(filename)) {
					return;
				}
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
				
				targetAccessionsByPeptide.put(match.getPeptideSeq(), match.getAccessions());
				
				String psmID=PercolatorPeptide.getPSMID(rawFile.getName(), retentionTime, Optional.empty(), false, peptideModSeq, precursorCharge);
				PercolatorPeptide pep=new PercolatorPeptide(psmID, PSMData.accessionsToString(accessionSet), 
						qValue, PEP, params.getAAConstants());
				scoredPeptides.add(pep);
			}
			
			@Override
			public void cleanup() {
			}
		});
		
		Logger.logLine("Found "+matches.size()+" peaks to integrate corresponding to "+targetAccessionsByPeptide.size()+" unique peptides");
		
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
		HashMap<SearchJobData, TObjectFloatHashMap<String>> rtByPeptideModSeq=new HashMap<SearchJobData, TObjectFloatHashMap<String>>();
		
		LibraryPeakLocationInferrer inferrer=new LibraryPeakLocationInferrer(selectedPredictions.values(), rtByPeptideModSeq, params);
		Logger.logLine("Finished predictions! Able to capture "+selectedPredictions.size()+" total predictions");
		
		StripeFileInterface stripeFile=StripeFileGenerator.getFile(rawFile, params, true);
		
		PeptideQuantExtractor extractor=new PeptideQuantExtractor(progress, stripeFile, params);
		ArrayList<IntegratedLibraryEntry> entries=extractor.extractPeptides(matches, Optional.of(inferrer), false);

		LibraryFile newLibrary=new LibraryFile();
		newLibrary.openFile();
		Logger.logLine("Writing library file "+newLibrary.getName());
		newLibrary.dropIndices();
		
		newLibrary.addTIC(stripeFile);
		
		newLibrary.addIntegratedEntries(!integratePrecursors, entries, Optional.empty(), Optional.empty(), params);		
		
		newLibrary.addProteinsFromEntries(targetAccessionsByPeptide, new HashMap<String, HashSet<String>>());
		newLibrary.addTargetDecoyPeptides(scoredPeptides, decoyPeptides);
		Pair<ArrayList<PercolatorProteinGroup>, ArrayList<PercolatorProteinGroup>> targetDecoyProteins=ParsimonyProteinGrouper.groupProteins(scoredPeptides, decoyPeptides, params.getPercolatorProteinThreshold(), params.getAAConstants());

		Logger.logLine("Writing global target/decoy proteins: "+targetDecoyProteins.x.size()+"/"+targetDecoyProteins.y.size());
		newLibrary.addTargetDecoyProteins(rawFile.getName(), targetDecoyProteins.x, targetDecoyProteins.y);

		newLibrary.createIndices();
		newLibrary.saveAsFile(newLibraryFile);
	}
	
	private static String getPeptideKey(String peptideModSeq, byte precursorCharge) {
		return peptideModSeq+"_"+precursorCharge;
	}
}
