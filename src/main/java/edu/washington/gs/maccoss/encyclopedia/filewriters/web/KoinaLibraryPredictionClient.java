package edu.washington.gs.maccoss.encyclopedia.filewriters.web;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.zip.DataFormatException;

import org.apache.commons.lang3.StringUtils;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
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

public class KoinaLibraryPredictionClient {
	private static final int BATCH_SIZE=1000;
	private final ArrayList<KoinaFeaturePredictionModel> models=new ArrayList<KoinaFeaturePredictionModel>();
	
	public static void main(String[] args) throws Exception {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		
		ArrayList<KoinaPrecursor> precursors=new ArrayList<KoinaPrecursor>();
		precursors.add(new KoinaPrecursor("LGGNEQVCR", 25f, (byte)2));
		precursors.add(new KoinaPrecursor("GAGSSEPVTGLDAK", 25f, (byte)2));
		
		AminoAcidConstants constants = new AminoAcidConstants();
		ArrayList<KoinaFeaturePredictionModel> models=new ArrayList<KoinaFeaturePredictionModel>();
		models.add(new Prosit2023timsTOFModel());
		models.add(new DeepLCHelaRTModel());
		models.add(new AlphaPeptDeepIMSModel());
		KoinaLibraryPredictionClient client=new KoinaLibraryPredictionClient(models);
		runKoinaOnBatch(client, precursors);
		
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
	
	public static void writeLibrary(ArrayList<KoinaFeaturePredictionModel> models, String libFileName, LibraryFile inputLibrary, int defaultNCE, byte defaultCharge, boolean adjustNCEForDIA, boolean addDecoys, SearchParameters params, ProgressIndicator progress) throws FileNotFoundException {
		libFileName = checkLibName(libFileName, inputLibrary.getFile(), null, defaultNCE, defaultCharge);
		try {
			LibraryFile library=new LibraryFile();
			library.openFile();
			Logger.logLine("Starting to build Prosit Library: "+libFileName);
			
			progress.update("Reading peptides from Library", 0.01f);
			ArrayList<KoinaPrecursor> allPeptides=getPeptidesFromLibrary(inputLibrary, defaultNCE, defaultCharge, adjustNCEForDIA, addDecoys, params.getAAConstants());
			
			SubProgressIndicator subProgress=new SubProgressIndicator(progress, 0.98f);
			int total=loadPredictionsIntoLibrary(allPeptides, library, models, subProgress);
			
			Logger.logLine("Processed "+total+" peptides, saving...");
			progress.update("Writing library to disk", 0.01f);
			library.saveAsFile(new File(libFileName));
			library.close();
			
			Logger.logLine("Finished writing "+total+" peptides to Prosit library!");
		} catch (IOException | SQLException | DataFormatException | InterruptedException e) {
			throw new EncyclopediaException("Unexpected error with Koina", e);
		}
	}
	
	public static void writeLibrary(ArrayList<KoinaFeaturePredictionModel> models, String libFileName, File fasta, DigestionEnzyme enzyme, int defaultNCE, byte defaultCharge, byte minCharge, byte maxCharge, int maxMissedCleavages, Range mzRange, boolean adjustNCEForDIA, boolean addDecoys, SearchParameters params, ProgressIndicator progress) throws FileNotFoundException {
		libFileName = checkLibName(libFileName, fasta, enzyme, defaultNCE, defaultCharge);
		try {
			LibraryFile library=new LibraryFile();
			library.openFile();
			Logger.logLine("Starting to build Prosit Library: "+libFileName);

			progress.update("Reading peptides from FASTA", 0.01f);
			ArrayList<KoinaPrecursor> allPeptides=getPeptidesFromFASTA(fasta, enzyme, minCharge, maxCharge, maxMissedCleavages, mzRange, defaultNCE, defaultCharge, adjustNCEForDIA, addDecoys, params.getAAConstants());

			SubProgressIndicator subProgress=new SubProgressIndicator(progress, 0.98f);
			int total=loadPredictionsIntoLibrary(allPeptides, library, models, subProgress);
			
			Logger.logLine("Processed "+total+" peptides, saving...");
			progress.update("Writing library to disk", 0.01f);
			library.saveAsFile(new File(libFileName));
			library.close();
			
			Logger.logLine("Finished writing "+total+" peptides to Prosit library!");
		} catch (IOException | SQLException | InterruptedException e) {
			throw new EncyclopediaException("Unexpected error with Koina", e);
		}
	}
	
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
	
	private static ArrayList<KoinaPrecursor> getPeptidesFromLibrary(LibraryFile inputLibrary, int defaultNCE, byte defaultCharge, boolean adjustNCEForDIA, boolean addDecoys, AminoAcidConstants constants) throws IOException, SQLException, DataFormatException {
		HashMap<KoinaPrecursor, KoinaPrecursor> allPeptides=new HashMap<>();

		AminoAcidConstants aminoAcidConstants = new AminoAcidConstants();
		ArrayList<LibraryEntry> allEntries = inputLibrary.getAllEntries(false, aminoAcidConstants);
		
		for (LibraryEntry entry : allEntries) {
			String seq=entry.getPeptideSeq();
			byte pepCharge=entry.getPrecursorCharge();
			
			KoinaPrecursor precursor=getKoinaPeptide(seq, pepCharge, defaultNCE, defaultCharge, adjustNCEForDIA, constants);
			
			if (precursor!=null) {
				KoinaPrecursor previous=allPeptides.get(precursor);
				if (previous==null) {
					precursor.addAccessions(entry.getAccessions());
					allPeptides.put(precursor, precursor);
				} else {
					previous.addAccessions(entry.getAccessions());
				}
				if (addDecoys) {
					String reverse=PeptideUtils.reverse(seq, constants);
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
	
	private static ArrayList<KoinaPrecursor> getPeptidesFromFASTA(File fasta, DigestionEnzyme enzyme, byte minCharge, byte maxCharge, int maxMissedCleavages, Range mzRange, int defaultNCE, byte defaultCharge, boolean adjustNCEForDIA, boolean addDecoys, AminoAcidConstants constants) {
		HashMap<KoinaPrecursor, KoinaPrecursor> allPeptides=new HashMap<>();
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(fasta, parameters);
		AminoAcidConstants aminoAcidConstants = new AminoAcidConstants();
		
		for (FastaEntryInterface entry : entries) {
			ArrayList<FastaPeptideEntry> peptidesInProtein=enzyme.digestProtein(entry, 7, 30, maxMissedCleavages, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
			for (FastaPeptideEntry pep : peptidesInProtein) {
				for (byte pepCharge = minCharge; pepCharge <=maxCharge; pepCharge++) {
					String seq=pep.getSequence();
					double pepMass=aminoAcidConstants.getMass(seq)+MassConstants.oh2;
					double pepChargedMass=(pepMass+MassConstants.protonMass*pepCharge)/pepCharge;

					if (mzRange.contains(pepChargedMass)) {
						if (seq.indexOf('B')>=0||seq.indexOf('J')>=0||seq.indexOf('O')>=0||seq.indexOf('U')>=0||seq.indexOf('X')>=0||seq.indexOf('Z')>=0||seq.indexOf('*')>=0) {
							continue;
						} else {
							KoinaPrecursor precursor=getKoinaPeptide(seq, pepCharge, defaultNCE, defaultCharge, adjustNCEForDIA, constants);
							
							if (precursor!=null) {
								KoinaPrecursor previous=allPeptides.get(precursor);
								if (previous==null) {
									precursor.addAccession(entry.getAccession());
									allPeptides.put(precursor, precursor);
								} else {
									previous.addAccession(entry.getAccession());
								}
								if (addDecoys) {
									String reverse=PeptideUtils.reverse(seq, constants);
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
			}
		}
		
		return new ArrayList<>(allPeptides.keySet());
	}

	private static KoinaPrecursor getKoinaPeptide(String seq, byte precursorCharge, int defaultNCE, byte defaultCharge, boolean adjustNCEForDIA, AminoAcidConstants constants) {
			// Prosit doesn't support charge >6
			if (precursorCharge<1&&precursorCharge>6) {
				return null;
			}
			
			// remove peptides that don't match PROSIT limitations:
			if (seq.indexOf('B')>=0||seq.indexOf('J')>=0||seq.indexOf('O')>=0||seq.indexOf('U')>=0||seq.indexOf('X')>=0||seq.indexOf('Z')>=0||seq.indexOf('*')>=0) {
				return null;
			}
			if (seq.length()<7||seq.length()>30) {
				return null;
			}
			
			float nce = adjustNCEForDIA?PrositCSVWriter.convertNCE(defaultNCE, precursorCharge, defaultCharge):defaultNCE;
			return new KoinaPrecursor(seq, nce, precursorCharge);
		
	}

	private static HashSet<KoinaPrecursor> getKoinaPeptides(int defaultNCE, byte defaultCharge, boolean adjustNCEForDIA, boolean addDecoys, HashSet<PeptidePrecursor> allPeptides, AminoAcidConstants constants) {
		HashSet<KoinaPrecursor> writablePeptides=new HashSet<>();
		
		for (PeptidePrecursor peptidePrecursor : allPeptides) {
			// Prosit only supports unmodified peptides:
			String seq = peptidePrecursor.getPeptideSeq();
			byte precursorCharge = peptidePrecursor.getPrecursorCharge();

			// Prosit doesn't support charge >6
			if (precursorCharge<1&&precursorCharge>6) {
				continue;
			}
			
			// remove peptides that don't match PROSIT limitations:
			if (seq.indexOf('B')>=0||seq.indexOf('J')>=0||seq.indexOf('O')>=0||seq.indexOf('U')>=0||seq.indexOf('X')>=0||seq.indexOf('Z')>=0||seq.indexOf('*')>=0) {
				continue;
			}
			if (seq.length()<7||seq.length()>30) {
				continue;
			}
			
			float nce = adjustNCEForDIA?PrositCSVWriter.convertNCE(defaultNCE, precursorCharge, defaultCharge):defaultNCE;
			writablePeptides.add(new KoinaPrecursor(seq, nce, precursorCharge));
			
			if (addDecoys) {
				String reverse=PeptideUtils.reverse(seq, constants);
				writablePeptides.add(new KoinaPrecursor(reverse, nce, precursorCharge));
			}
		}
		return writablePeptides;
	}

	private KoinaLibraryPredictionClient(ArrayList<KoinaFeaturePredictionModel> models) {
		this.models.addAll(models);
	}
	
	private static int loadPredictionsIntoLibrary(ArrayList<KoinaPrecursor> peptides, LibraryFile library, ArrayList<KoinaFeaturePredictionModel> models, ProgressIndicator progress) throws IOException, SQLException, InterruptedException {
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
				runKoinaOnBatch(client, subList);

				for (KoinaPrecursor precursor : subList) {
					count++;
					AnnotatedLibraryEntry entry=precursor.toEntry(constants, params);
					returnList.add(entry);
				}
				
			} catch (Exception e) {
				try {
					Logger.errorLine("Ran into a Koina error, trying a second time!");
					runKoinaOnBatch(client, subList);

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

	private static void runKoinaOnBatch(KoinaLibraryPredictionClient client, List<KoinaPrecursor> subList)
			throws InterruptedException {
		try {
			ArrayList<Thread> threads=new ArrayList<Thread>();
			for (KoinaFeaturePredictionModel model : client.models) {
				Thread t=new Thread(new Runnable() {
					@Override
					public void run() {
						model.updatePeptides(subList);
					}
				});
				t.start();
				threads.add(t);
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
