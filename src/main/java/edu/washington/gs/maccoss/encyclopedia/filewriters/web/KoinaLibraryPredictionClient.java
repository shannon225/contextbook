package edu.washington.gs.maccoss.encyclopedia.filewriters.web;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Optional;
import java.util.zip.DataFormatException;

import org.apache.commons.lang3.StringUtils;

import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;
import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.NonstandardAminoAcidException;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PrositCSVWriter;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation.Prosit2020HCDModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation.Prosit2023timsTOFModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.ims.AlphaPeptDeepIMSModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.ims.IM2DeepIMSModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt.DeepLCHelaRTModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt.Prosit2019RTModel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.RunnableWithExceptions;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;
import gnu.trove.map.hash.TCharDoubleHashMap;

/**
 * KoinaLibraryPredictionClient is responsible for managing the execution of feature prediction
 * models and building a Koina library from provided precursors.
 */
public class KoinaLibraryPredictionClient {
	public static final String HTTPS_KOINA_WILHELMLAB_ORG_443 = "https://koina.wilhelmlab.org:443/";
	
	private static final int BATCH_SIZE=5000;
	private static final int MICRO_BATCH_SIZE=1000;
	private final ArrayList<KoinaFeaturePredictionModel> models=new ArrayList<KoinaFeaturePredictionModel>();
	
	/**
     * Main method to demonstrate how to create a KoinaLibraryPredictionClient and predict features
     * for a batch of precursors.
     *
     * @param args command-line arguments (not used).
     * @throws Exception if any error occurs during execution.
     */
	public static void main(String[] args) throws Exception {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		
		ArrayList<KoinaPrecursor> precursors=new ArrayList<KoinaPrecursor>();
		precursors.add(new KoinaPrecursor(AminoAcidEncoding.getAAs("LGGNEQVCR", params.getAAConstants()), 25f, (byte)2));
		precursors.add(new KoinaPrecursor(AminoAcidEncoding.getAAs("GAGSSEPVTGLDAK", params.getAAConstants()), 25f, (byte)2));
		
		AminoAcidConstants constants = new AminoAcidConstants();
		ArrayList<KoinaFeaturePredictionModel> models=new ArrayList<KoinaFeaturePredictionModel>();
		models.add(new Prosit2023timsTOFModel());
		models.add(new DeepLCHelaRTModel());
		models.add(new AlphaPeptDeepIMSModel());
		KoinaLibraryPredictionClient client=new KoinaLibraryPredictionClient(models);
		runKoinaOnBatch(client, precursors, KoinaLibraryPredictionClient.HTTPS_KOINA_WILHELMLAB_ORG_443);
		
		for (KoinaPrecursor precursor : precursors) {
			AnnotatedLibraryEntry entry=precursor.toEntry(constants, params);
			Charter.launchChart(entry);
		}
	}

	public static ArrayList<KoinaFeaturePredictionModel> getDefaultModels() {
		ArrayList<KoinaFeaturePredictionModel> models=new ArrayList<KoinaFeaturePredictionModel>();
		models.add(new Prosit2020HCDModel());
		models.add(new Prosit2019RTModel());
		models.add(new IM2DeepIMSModel());
		return models;
	}
	
	/**
     * Writes a Koina library using the provided feature prediction models and input library.
     *
     * @param models        the list of feature prediction models to use.
     * @param libFileName   the name of the output library file.
     * @param inputLibrary  the input library file containing peptides.
     * @param defaultNCE    the default normalized collision energy (NCE) to apply.
     * @param defaultCharge the default charge state to use for precursors.
     * @param adjustNCEForDIA whether to adjust NCE for DIA experiments.
     * @param addDecoys     whether to add decoy entries.
     * @param params        the search parameters object.
     * @param progress      a progress indicator to track the library building process.
     * @throws FileNotFoundException if the input library file is not found.
     */
	public static void writeLibrary(String baseURL, ArrayList<KoinaFeaturePredictionModel> models, String libFileName, LibraryFile inputLibrary, int defaultNCE, byte defaultCharge, boolean adjustNCEForDIA, boolean addDecoys, SearchParameters params, ProgressIndicator progress) throws FileNotFoundException {
		libFileName = checkLibName(libFileName, inputLibrary.getFile(), null, defaultNCE, defaultCharge);
		try {
			LibraryFile library=new LibraryFile();
			library.openFile();
			Logger.logLine("Starting to build Koina Library: "+libFileName+", models used:");
			for (KoinaFeaturePredictionModel model : models) {
				Logger.logLine("\t"+model.getName()+" ("+model.getModelType()+")");
			}
			Logger.logLine("Using Server: ["+baseURL+"]");
			
			progress.update("Reading peptides from Library", 0.01f);
			ArrayList<KoinaPrecursor> allPeptides=getPeptidesFromLibrary(models, inputLibrary, defaultNCE, defaultCharge, adjustNCEForDIA, addDecoys, params.getAAConstants());
			
			SubProgressIndicator subProgress=new SubProgressIndicator(progress, 0.98f);
			int total=loadPredictionsIntoLibrary(baseURL, allPeptides, library, models, subProgress);
			
			Logger.logLine("Processed "+total+" peptides, saving...");
			progress.update("Writing library to disk", 0.01f);
			library.saveAsFile(new File(libFileName));
			library.close();
			
			Logger.logLine("Finished writing "+total+" peptides to Koina library!");
		} catch (IOException | SQLException | DataFormatException | InterruptedException e) {
			throw new EncyclopediaException("Unexpected error with Koina", e);
		}
	}
	
	/**
     * Writes a Koina library using the provided feature prediction models and peptides digested
     * from a FASTA file.
     *
     * @param models        the list of feature prediction models to use.
     * @param libFileName   the name of the output library file.
     * @param fasta         the FASTA file containing protein sequences.
     * @param enzyme        the digestion enzyme used for peptide digestion.
     * @param defaultNCE    the default normalized collision energy (NCE) to apply.
     * @param defaultCharge the default charge state to use for precursors.
     * @param minCharge     the minimum charge state for peptides.
     * @param maxCharge     the maximum charge state for peptides.
     * @param maxMissedCleavages the maximum number of missed cleavages allowed.
     * @param mzRange       the mass-to-charge ratio range for valid peptides.
     * @param adjustNCEForDIA whether to adjust NCE for DIA experiments.
     * @param addDecoys     whether to add decoy entries.
     * @param params        the search parameters object.
     * @param progress      a progress indicator to track the library building process.
     * @throws FileNotFoundException if the FASTA file is not found.
     */
	public static void writeLibrary(String baseURL, ArrayList<KoinaFeaturePredictionModel> models, String libFileName, File fasta, DigestionEnzyme enzyme, int defaultNCE, byte defaultCharge, byte minCharge, byte maxCharge, int maxMissedCleavages, Range mzRange, boolean adjustNCEForDIA, boolean addDecoys, SearchParameters params, ProgressIndicator progress) throws FileNotFoundException {
		libFileName = checkLibName(libFileName, fasta, enzyme, defaultNCE, defaultCharge);
		try {
			LibraryFile library=new LibraryFile();
			library.openFile();
			Logger.logLine("Starting to build Koina Library: "+libFileName+", models used:");
			for (KoinaFeaturePredictionModel model : models) {
				Logger.logLine("\t"+model.getName()+" ("+model.getModelType()+")");
			}
			Logger.logLine("Using Server: ["+baseURL+"]");

			progress.update("Reading peptides from FASTA", 0.01f);
			ArrayList<KoinaPrecursor> allPeptides=getPeptidesFromFASTA(models, fasta, enzyme, minCharge, maxCharge, maxMissedCleavages, mzRange, defaultNCE, defaultCharge, adjustNCEForDIA, addDecoys, params.getAAConstants());

			SubProgressIndicator subProgress=new SubProgressIndicator(progress, 0.98f);
			int total=loadPredictionsIntoLibrary(baseURL, allPeptides, library, models, subProgress);
			
			Logger.logLine("Processed "+total+" peptides, saving...");
			progress.update("Writing library to disk", 0.01f);
			library.saveAsFile(new File(libFileName));
			library.close();
			
			Logger.logLine("Finished writing "+total+" peptides to Koina library!");
		} catch (IOException | SQLException | InterruptedException e) {
			throw new EncyclopediaException("Unexpected error with Koina", e);
		}
	}
	
	/**
     * Constructs a valid library filename if none is provided or validates the provided name.
     */
	private static String checkLibName(String csvFileName, File fasta, DigestionEnzyme enzyme, int defaultNCE, byte defaultCharge) {
		String fileName;
		if (null==csvFileName||StringUtils.isBlank(csvFileName)) {
			String enzymeText=enzyme==null?"":("."+enzyme.getPercolatorName());
			fileName = fasta.getAbsolutePath() + enzymeText + ".z" + defaultCharge + "_nce" + defaultNCE + LibraryFile.DLIB;
		} else {
			fileName = csvFileName;
		}
		return fileName;
	}
	
	/**
     * Extracts valid peptides from a library file, ensuring they are compatible with all models.
     */
	private static ArrayList<KoinaPrecursor> getPeptidesFromLibrary(ArrayList<KoinaFeaturePredictionModel> models, LibraryFile inputLibrary, int defaultNCE, byte defaultCharge, boolean adjustNCEForDIA, boolean addDecoys, AminoAcidConstants constants) throws IOException, SQLException, DataFormatException {
		HashMap<KoinaPrecursor, KoinaPrecursor> allPeptides=new HashMap<>();

		AminoAcidConstants aminoAcidConstants = new AminoAcidConstants();
		ArrayList<LibraryEntry> allEntries = inputLibrary.getAllEntries(false, aminoAcidConstants);
		
		for (LibraryEntry entry : allEntries) {
			byte pepCharge=entry.getPrecursorCharge();

			boolean passes=true;
			AminoAcidEncoding[] aas=null;
			try {
				aas=AminoAcidEncoding.getAAs(entry.getPeptideModSeq(), constants);
				for (KoinaFeaturePredictionModel model : models) {
					if (!model.canModelPeptide(aas, pepCharge)) {
						passes=false;
						break;
					}
				}
			} catch (NonstandardAminoAcidException e) {
				passes=false;
			}
			if (!passes) continue;
			
			KoinaPrecursor precursor=getKoinaPeptide(aas, pepCharge, defaultNCE, defaultCharge, adjustNCEForDIA, constants);
			
			if (precursor!=null) {
				KoinaPrecursor previous=allPeptides.get(precursor);
				if (previous==null) {
					precursor.addAccessions(entry.getAccessions());
					allPeptides.put(precursor, precursor);
				} else {
					previous.addAccessions(entry.getAccessions());
				}
				if (addDecoys) {
					AminoAcidEncoding[] reverse=AminoAcidEncoding.reverse(aas.clone());
					KoinaPrecursor revPrecursor = new KoinaPrecursor(reverse, precursor.getNCE(), precursor.getCharge());
					
					previous=allPeptides.get(revPrecursor);
					if (previous==null) {
						for (String accession : entry.getAccessions()) {
							revPrecursor.addAccession(LibraryEntry.DECOY_STRING+accession);
						}
						allPeptides.put(revPrecursor, revPrecursor);
					} else {
						for (String accession : entry.getAccessions()) {
							previous.addAccession(LibraryEntry.DECOY_STRING+accession);
						}
					}
				}
			}
		}
		
		return new ArrayList<>(allPeptides.keySet());
	}

    /**
     * Extracts valid peptides from a FASTA file, ensuring they are compatible with all models.
     */
	private static ArrayList<KoinaPrecursor> getPeptidesFromFASTA(ArrayList<KoinaFeaturePredictionModel> models, File fasta, DigestionEnzyme enzyme, byte minCharge, byte maxCharge, int maxMissedCleavages, Range mzRange, int defaultNCE, byte defaultCharge, boolean adjustNCEForDIA, boolean addDecoys, AminoAcidConstants constants) {
		HashMap<KoinaPrecursor, KoinaPrecursor> allPeptides=new HashMap<>();
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(fasta, parameters);
		
		// assumes C+57
		AminoAcidConstants aminoAcidConstants = new AminoAcidConstants();
		
		for (FastaEntryInterface entry : entries) {
			ArrayList<FastaPeptideEntry> peptidesInProtein=enzyme.digestProtein(entry, 7, 30, maxMissedCleavages, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
			for (FastaPeptideEntry pep : peptidesInProtein) {
				for (byte pepCharge = minCharge; pepCharge <=maxCharge; pepCharge++) {
					String peptideModSeq=aminoAcidConstants.toPeptideModSeq(pep.getSequence());
					
					double pepMass=aminoAcidConstants.getMass(peptideModSeq)+MassConstants.oh2;
					double pepChargedMass=(pepMass+MassConstants.protonMass*pepCharge)/pepCharge;

					if (!mzRange.contains(pepChargedMass)) {
						continue;
					}

					boolean passes=true;
					AminoAcidEncoding[] aas=null;
					try {
						aas=AminoAcidEncoding.getAAs(peptideModSeq, constants);
						for (KoinaFeaturePredictionModel model : models) {
							if (!model.canModelPeptide(aas, pepCharge)) {
								passes=false;
								break;
							}
						}
					} catch (NonstandardAminoAcidException e) {
						passes=false;
					}
					if (!passes) continue;

					KoinaPrecursor precursor=getKoinaPeptide(aas, pepCharge, defaultNCE, defaultCharge, adjustNCEForDIA, constants);
					
					if (precursor!=null) {
						KoinaPrecursor previous=allPeptides.get(precursor);
						if (previous==null) {
							precursor.addAccession(entry.getAccession());
							allPeptides.put(precursor, precursor);
						} else {
							previous.addAccession(entry.getAccession());
						}
						if (addDecoys) {
							AminoAcidEncoding[] reverse=AminoAcidEncoding.reverse(aas.clone());
							KoinaPrecursor revPrecursor = new KoinaPrecursor(reverse, precursor.getNCE(), precursor.getCharge());
						
							previous=allPeptides.get(revPrecursor);
							if (previous==null) {
								revPrecursor.addAccession(LibraryEntry.DECOY_STRING+entry.getAccession());
								allPeptides.put(revPrecursor, revPrecursor);
							} else {
								previous.addAccession(LibraryEntry.DECOY_STRING+entry.getAccession());
							}
						}
					}
				}
			}
		}
		
		return new ArrayList<>(allPeptides.keySet());
	}

    /**
     * Creates a KoinaPrecursor from an amino acid sequence.
     */
	private static KoinaPrecursor getKoinaPeptide(AminoAcidEncoding[] aas, byte precursorCharge, int defaultNCE, byte defaultCharge, boolean adjustNCEForDIA, AminoAcidConstants constants) {
			float nce = adjustNCEForDIA?PrositCSVWriter.convertNCE(defaultNCE, precursorCharge, defaultCharge):defaultNCE;
			return new KoinaPrecursor(aas, nce, precursorCharge);
		
	}
	
	/**
	 * Constructs a new KoinaLibraryPredictionClient with the specified list of feature prediction models.
	 */
	private KoinaLibraryPredictionClient(ArrayList<KoinaFeaturePredictionModel> models) {
		this.models.addAll(models);
	}
	
	/**
	 * Loads the predictions from the feature prediction models into the library by processing the provided peptides.
	 */
	private static int loadPredictionsIntoLibrary(String baseURL, ArrayList<KoinaPrecursor> peptides, LibraryFile library, ArrayList<KoinaFeaturePredictionModel> models, ProgressIndicator progress) throws IOException, SQLException, InterruptedException {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		AminoAcidConstants constants = new AminoAcidConstants();
		KoinaLibraryPredictionClient client=new KoinaLibraryPredictionClient(models);
		
		library.dropIndices();
		
		int count=0;
		int start=0;
		// track existing thread and runnable for finishing and  exceptions
		Thread prevThread=null;
		RunnableWithExceptions prevRunnable=null;
		while (true) {
			if (start>=peptides.size()) {
				break;
			}
			
			int stop=Math.min(peptides.size(), start+BATCH_SIZE);
			List<KoinaPrecursor> subList=peptides.subList(start, stop);
			
			// try koina twice in a row before failing
			ArrayList<LibraryEntry> returnList=new ArrayList<LibraryEntry>();
			try {
				runKoinaOnBatch(client, subList, baseURL);

				for (KoinaPrecursor precursor : subList) {
					count++;
					AnnotatedLibraryEntry entry=precursor.toEntry(constants, params);
					returnList.add(entry);
				}
				
			} catch (Exception e) {
				try {
					Logger.errorLine("Ran into a Koina error, trying a second time!");
					runKoinaOnBatch(client, subList, baseURL);

					returnList.clear();
					for (KoinaPrecursor precursor : subList) {
						count++;
						AnnotatedLibraryEntry entry=precursor.toEntry(constants, params);
						returnList.add(entry);
					}
					
				} catch (Exception e2) {
					Logger.errorLine("Ran into a second Koina error, failing!");
					throw new EncyclopediaException(e2);
				}
			}

			// Join previous because we can't have two of these threads running at the same time!
			// This is because they write to the library, which is not thread safe
			if (prevThread!=null) {
				prevThread.join();
				Optional<Exception> maybeException = prevRunnable.getException();
				if (maybeException.isPresent()) {
					throw new EncyclopediaException(maybeException.get());
				}
			}
			
			prevRunnable=new RunnableWithExceptions() {
				@Override
				public void run() {
					try {
						library.addEntries(returnList);
						library.addProteinsFromEntries(returnList);
					} catch (IOException | SQLException e) {
						setException(e);
					}
				}
				
			};
			prevThread=new Thread(prevRunnable);
			prevThread.start();
			
			Logger.logLine("Processed "+count+" of "+peptides.size());
			start=stop;
			progress.update("Processed "+count+" of "+peptides.size(), count/(float)peptides.size());
		}
		
		// finish previous writing thread if still dangling
		if (prevThread!=null) {
			prevThread.join();
			Optional<Exception> maybeException = prevRunnable.getException();
			if (maybeException.isPresent()) {
				throw new EncyclopediaException(maybeException.get());
			}
		}
		
		library.createIndices();
		return count;
	}


	/**
	 * Runs the Koina feature prediction models on a batch of KoinaPrecursor peptides.
	 * This method creates and starts threads for each model to update the peptides in parallel.
	 */
	private static void runKoinaOnBatch(KoinaLibraryPredictionClient client, List<KoinaPrecursor> peptides, String baseURL)
			throws InterruptedException {
		try {
			ArrayList<Thread> threads=new ArrayList<Thread>();
			int start=0;
			while (true) {
				if (start>=peptides.size()) {
					break;
				}
				
				int stop=Math.min(peptides.size(), start+MICRO_BATCH_SIZE);
				List<KoinaPrecursor> subList=peptides.subList(start, stop);
			
			
				for (KoinaFeaturePredictionModel model : client.models) {
					Thread t=new Thread(new Runnable() {
						@Override
						public void run() {
							model.updatePeptides(subList, baseURL);
						}
					});
					t.start();
					threads.add(t);
				}
				
				start=stop;
			}

			for (Thread t : threads) {
				t.join(60*1000);
			}

			boolean allSucceeded=true;
			for (Thread t : threads) {
				if (t.isAlive()) {
					allSucceeded=false;
					break;
				}
			}

			if (!allSucceeded) {
				throw new EncyclopediaException("Timed out on Koina job");
			}

		} catch (Exception e) { 
			throw new EncyclopediaException("Failed query on Koina Job", e);
		}
	}
}
