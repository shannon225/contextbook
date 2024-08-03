package edu.washington.gs.maccoss.encyclopedia.filewriters.web;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.zip.DataFormatException;

import org.apache.commons.lang3.StringUtils;
import org.json.JSONArray;
import org.json.JSONObject;

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
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.MatrixMath;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class KoinaLibraryPredictionClient {
	private static final int BATCH_SIZE=512;
	boolean isHCD;
	
	public static void main(String[] args) {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		
		ArrayList<KoinaPrecursor> precursors=new ArrayList<KoinaPrecursor>();
		precursors.add(new KoinaPrecursor("LGGNEQVCR", 25f, (byte)2));
		precursors.add(new KoinaPrecursor("GAGSSEPVTGLDAK", 25f, (byte)2));
		
		AminoAcidConstants constants = new AminoAcidConstants();
		KoinaLibraryPredictionClient client=new KoinaLibraryPredictionClient(true);
		
		client.postPeptideIntensities(precursors);
		client.postPeptideRTs(precursors);
		
		for (KoinaPrecursor precursor : precursors) {
			AnnotatedLibraryEntry entry=precursor.toEntry(constants, params);
			Charter.launchChart(entry);
		}
	}
	
	public static void writeLibrary(String libFileName, LibraryFile inputLibrary, int defaultNCE, byte defaultCharge, boolean adjustNCEForDIA, boolean addDecoys, SearchParameters params, ProgressIndicator progress) throws FileNotFoundException {
		libFileName = checkLibName(libFileName, inputLibrary.getFile(), null, defaultNCE, defaultCharge);
		try {
			LibraryFile library=new LibraryFile();
			library.openFile();
			Logger.logLine("Starting to build Prosit Library: "+libFileName);
			
			progress.update("Reading peptides from Library", 0.01f);
			ArrayList<KoinaPrecursor> allPeptides=getPeptidesFromLibrary(inputLibrary, defaultNCE, defaultCharge, adjustNCEForDIA, addDecoys, params.getAAConstants());
			
			SubProgressIndicator subProgress=new SubProgressIndicator(progress, 0.98f);
			int total=loadPredictionsIntoLibrary(allPeptides, library, subProgress);
			
			Logger.logLine("Processed "+total+" peptides, saving...");
			progress.update("Writing library to disk", 0.01f);
			library.saveAsFile(new File(libFileName));
			library.close();
			
			Logger.logLine("Finished writing "+total+" peptides to Prosit library!");
		} catch (IOException | SQLException | DataFormatException e) {
			
		}
	}
	
	public static void writeLibrary(String libFileName, File fasta, DigestionEnzyme enzyme, int defaultNCE, byte defaultCharge, byte minCharge, byte maxCharge, int maxMissedCleavages, Range mzRange, boolean adjustNCEForDIA, boolean addDecoys, SearchParameters params, ProgressIndicator progress) throws FileNotFoundException {
		libFileName = checkLibName(libFileName, fasta, enzyme, defaultNCE, defaultCharge);
		try {
			LibraryFile library=new LibraryFile();
			library.openFile();
			Logger.logLine("Starting to build Prosit Library: "+libFileName);

			progress.update("Reading peptides from FASTA", 0.01f);
			ArrayList<KoinaPrecursor> allPeptides=getPeptidesFromFASTA(fasta, enzyme, minCharge, maxCharge, maxMissedCleavages, mzRange, defaultNCE, defaultCharge, adjustNCEForDIA, addDecoys, params.getAAConstants());

			SubProgressIndicator subProgress=new SubProgressIndicator(progress, 0.98f);
			int total=loadPredictionsIntoLibrary(allPeptides, library, subProgress);
			
			Logger.logLine("Processed "+total+" peptides, saving...");
			progress.update("Writing library to disk", 0.01f);
			library.saveAsFile(new File(libFileName));
			library.close();
			
			Logger.logLine("Finished writing "+total+" peptides to Prosit library!");
		} catch (IOException | SQLException e) {
			
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
					KoinaPrecursor revPrecursor = new KoinaPrecursor(reverse, precursor.nce, precursor.charge);
					
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
									KoinaPrecursor revPrecursor = new KoinaPrecursor(reverse, precursor.nce, precursor.charge);
									
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

	private KoinaLibraryPredictionClient(boolean isHCD) {
		this.isHCD=isHCD;
	}
	
	private static int loadPredictionsIntoLibrary(ArrayList<KoinaPrecursor> peptides, LibraryFile library, ProgressIndicator progress) throws IOException, SQLException {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		AminoAcidConstants constants = new AminoAcidConstants();
		KoinaLibraryPredictionClient client=new KoinaLibraryPredictionClient(true);
		
		library.dropIndices();
		
		int count=0;
		int start=0;
		while (true) {
			if (start>=peptides.size()) {
				break;
			}
			
			int stop=Math.min(peptides.size(), start+BATCH_SIZE);
			List<KoinaPrecursor> subList=peptides.subList(start, stop);
			
			try {
				runKoinaOnBatch(client, subList);
			} catch (Exception e) {

				try {
					Logger.errorLine("Ran into a Koina error, trying a second time!");
					runKoinaOnBatch(client, subList);
				} catch (Exception e2) {
					Logger.errorLine("Ran into a second Koina error, failing!");
					throw new EncyclopediaException(e2);
				}
			}

			ArrayList<LibraryEntry> returnList=new ArrayList<LibraryEntry>();
			for (KoinaPrecursor precursor : subList) {
				count++;
				AnnotatedLibraryEntry entry=precursor.toEntry(constants, params);
				returnList.add(entry);
			}
			library.addEntries(returnList);
			library.addProteinsFromEntries(returnList);
			
			Logger.logLine("Processed "+count+" of "+peptides.size());
			start=stop;
			progress.update("Processed "+count+" of "+peptides.size(), count/(float)peptides.size());
		}
		
		library.createIndices();
		return count;
	}

	private static void runKoinaOnBatch(KoinaLibraryPredictionClient client, List<KoinaPrecursor> subList)
			throws InterruptedException {
		try {
			Thread t1=new Thread(new Runnable() {
				@Override
				public void run() {
					client.postPeptideIntensities(subList);
				}
			});
			t1.start();
	
			Thread t2=new Thread(new Runnable() {
				@Override
				public void run() {
					client.postPeptideRTs(subList);
				}
			});
			t2.start();
	
			Thread t3=new Thread(new Runnable() {
				@Override
				public void run() {
					client.postPeptideIonMobilities(subList);
				}
			});
			t3.start();
			
			t1.join(60*1000);
			t2.join(60*1000);
			t3.join(60*1000);

			boolean allSucceeded=!t1.isAlive()&&!t2.isAlive()&&!t1.isAlive();
			
			if (!allSucceeded) {
				throw new EncyclopediaException("Timed out on Koina job");
			}

		} catch (Exception e) { 
			throw new EncyclopediaException("Failed on Koina Job", e);
		}
	}
	
	private URL getFragmentationURL() {
		try {
			if (isHCD) { 
				return new URL("https://koina.wilhelmlab.org/v2/models/Prosit_2020_intensity_HCD/infer");
			} else {
				return new URL("https://koina.wilhelmlab.org/v2/models/Prosit_2020_intensity_CID/infer");
			}
		} catch (MalformedURLException e) {
			throw new EncyclopediaException("Error getting Koina URL", e);
		}
	}
	
	private URL getIMSURL() {
		try {
			return new URL("https://koina.wilhelmlab.org:443/v2/models/IM2Deep/infer");
		} catch (MalformedURLException e) {
			throw new EncyclopediaException("Error getting Koina URL", e);
		}
	}
	
	private URL getRTURL() {
		try {
			return new URL("https://koina.wilhelmlab.org/v2/models/Prosit_2019_irt/infer");
		} catch (MalformedURLException e) {
			throw new EncyclopediaException("Error getting Koina URL", e);
		}
	}

	/**
	 * 
	 * @param peptides
	 * @return 
	 */
	private void postPeptideIonMobilities(List<KoinaPrecursor> peptides) {
		ArrayList<String> pepseqs=new ArrayList<String>();
		TFloatArrayList NCEs=new TFloatArrayList();
		TIntArrayList charges=new TIntArrayList();
		for (KoinaPrecursor pep : peptides) {
			pepseqs.add(pep.prositSequence);
			NCEs.add(pep.nce);
			charges.add(pep.charge);
		}
		
		try {
			URL url=getIMSURL();
			
			HttpURLConnection conn = (HttpURLConnection) url.openConnection();
	
	        // Set request method to POST
	        conn.setRequestMethod("POST");
	        conn.setRequestProperty("Content-Type", "application/json; utf-8");
	        conn.setRequestProperty("Accept", "application/json");
	        conn.setDoOutput(true);
	
	        // Create JSON body
	        JSONObject jsonBody = new JSONObject();
	        jsonBody.put("id", "EncyclopeDIA_query");
	
	        JSONArray inputsArray = new JSONArray();
	
	        JSONObject peptideSequences = new JSONObject();
	        peptideSequences.put("name", "peptide_sequences");
	        peptideSequences.put("shape", new JSONArray("["+peptides.size()+",1]"));
	        peptideSequences.put("datatype", "BYTES");
	        peptideSequences.put("data", new JSONArray(pepseqs.toArray(new String[0])));
	        inputsArray.put(peptideSequences);
	        
	        jsonBody.put("inputs", inputsArray);
	        JSONObject precursorCharges = new JSONObject();
	        precursorCharges.put("name", "precursor_charges");
	        precursorCharges.put("shape", new JSONArray("["+peptides.size()+",1]"));
	        precursorCharges.put("datatype", "INT32");
	        precursorCharges.put("data", new JSONArray(charges.toArray()));
	        inputsArray.put(precursorCharges);
	    	
	        // Write JSON body to request
	        try (OutputStream os = conn.getOutputStream()) {
	            byte[] input = jsonBody.toString().getBytes("utf-8");
	            os.write(input, 0, input.length);
	        }

            StringBuilder response = new StringBuilder();
	        // Read the response
	        try (BufferedReader br = new BufferedReader(
	                new InputStreamReader(conn.getInputStream(), "utf-8"))) {
	            String responseLine;
	            while ((responseLine = br.readLine()) != null) {
	                response.append(responseLine.trim());
	            }
	        }

	        parseIMSData(response.toString(), peptides);
	        
		} catch (IOException ioe) {
			throw new EncyclopediaException("IO error getting Koina result", ioe);
		}
	}

	/**
	 * 
	 * @param peptides
	 * @return 
	 */
	private void postPeptideRTs(List<KoinaPrecursor> peptides) {
		ArrayList<String> pepseqs=new ArrayList<String>();
		TFloatArrayList NCEs=new TFloatArrayList();
		TIntArrayList charges=new TIntArrayList();
		for (KoinaPrecursor pep : peptides) {
			pepseqs.add(pep.prositSequence);
			NCEs.add(pep.nce);
			charges.add(pep.charge);
		}
		
		try {
			URL url=getRTURL();
			
			HttpURLConnection conn = (HttpURLConnection) url.openConnection();
	
	        // Set request method to POST
	        conn.setRequestMethod("POST");
	        conn.setRequestProperty("Content-Type", "application/json; utf-8");
	        conn.setRequestProperty("Accept", "application/json");
	        conn.setDoOutput(true);
	
	        // Create JSON body
	        JSONObject jsonBody = new JSONObject();
	        jsonBody.put("id", "EncyclopeDIA_query");
	
	        JSONArray inputsArray = new JSONArray();
	
	        JSONObject peptideSequences = new JSONObject();
	        peptideSequences.put("name", "peptide_sequences");
	        peptideSequences.put("shape", new JSONArray("["+peptides.size()+",1]"));
	        peptideSequences.put("datatype", "BYTES");
	        peptideSequences.put("data", new JSONArray(pepseqs.toArray(new String[0])));
	        inputsArray.put(peptideSequences);
	        
	        jsonBody.put("inputs", inputsArray);
	    	
	        // Write JSON body to request
	        try (OutputStream os = conn.getOutputStream()) {
	            byte[] input = jsonBody.toString().getBytes("utf-8");
	            os.write(input, 0, input.length);
	        }

            StringBuilder response = new StringBuilder();
	        // Read the response
	        try (BufferedReader br = new BufferedReader(
	                new InputStreamReader(conn.getInputStream(), "utf-8"))) {
	            String responseLine;
	            while ((responseLine = br.readLine()) != null) {
	                response.append(responseLine.trim());
	            }
	        }

	        parseRTData(response.toString(), peptides);
	        
		} catch (IOException ioe) {
			throw new EncyclopediaException("IO error getting Koina result", ioe);
		}
	}
	
	/**
	 * 
	 * @param peptides
	 * @return 
	 */
	private void postPeptideIntensities(List<KoinaPrecursor> peptides) {
		ArrayList<String> pepseqs=new ArrayList<String>();
		TFloatArrayList NCEs=new TFloatArrayList();
		TIntArrayList charges=new TIntArrayList();
		for (KoinaPrecursor pep : peptides) {
			pepseqs.add(pep.prositSequence);
			NCEs.add(pep.nce);
			charges.add(pep.charge);
		}
		
		try {
			URL url=getFragmentationURL();
			
			HttpURLConnection conn = (HttpURLConnection) url.openConnection();
	
	        // Set request method to POST
	        conn.setRequestMethod("POST");
	        conn.setRequestProperty("Content-Type", "application/json; utf-8");
	        conn.setRequestProperty("Accept", "application/json");
	        conn.setDoOutput(true);
	
	        // Create JSON body
	        JSONObject jsonBody = new JSONObject();
	        jsonBody.put("id", "EncyclopeDIA_query");
	
	        JSONArray inputsArray = new JSONArray();
	
	        JSONObject peptideSequences = new JSONObject();
	        peptideSequences.put("name", "peptide_sequences");
	        peptideSequences.put("shape", new JSONArray("["+peptides.size()+",1]"));
	        peptideSequences.put("datatype", "BYTES");
	        peptideSequences.put("data", new JSONArray(pepseqs.toArray(new String[0])));
	
	        JSONObject precursorCharges = new JSONObject();
	        precursorCharges.put("name", "precursor_charges");
	        precursorCharges.put("shape", new JSONArray("["+peptides.size()+",1]"));
	        precursorCharges.put("datatype", "INT32");
	        precursorCharges.put("data", new JSONArray(charges.toArray()));
	
	        inputsArray.put(peptideSequences);
	        
	        if (isHCD) {
		        JSONObject collisionEnergies = new JSONObject();
		        collisionEnergies.put("name", "collision_energies");
		        collisionEnergies.put("shape", new JSONArray("["+peptides.size()+",1]"));
		        collisionEnergies.put("datatype", "FP32");
		        collisionEnergies.put("data", new JSONArray(NCEs.toArray()));
	        	inputsArray.put(collisionEnergies);
	        }
	        
	        inputsArray.put(precursorCharges);
	
	        jsonBody.put("inputs", inputsArray);
	
	        // Write JSON body to request
	        try (OutputStream os = conn.getOutputStream()) {
	            byte[] input = jsonBody.toString().getBytes("utf-8");
	            os.write(input, 0, input.length);
	        }

            StringBuilder response = new StringBuilder();
	        // Read the response
	        try (BufferedReader br = new BufferedReader(
	                new InputStreamReader(conn.getInputStream(), "utf-8"))) {
	            String responseLine;
	            while ((responseLine = br.readLine()) != null) {
	                response.append(responseLine.trim());
	            }
	        }
	
	        // Handle response code
	        int responseCode = conn.getResponseCode();
	        if (responseCode != HttpURLConnection.HTTP_OK) {
	            Logger.errorLine("Failed response code from Koina: "+responseCode);
	        }

	        parseFragmentationData(response.toString(), peptides);
	        
		} catch (IOException ioe) {
			throw new EncyclopediaException("IO error getting Koina result", ioe);
		}
	}

    private static void parseIMSData(String json, List<KoinaPrecursor> peptides) {
        json=removeGroupingBrackets(json);
        String[] pairs = splitJsonElements(json);
        
        String outputs=getSelectedElement(pairs, "\"outputs\":");
        String[] keyValue = outputs.split(":", 2);
        String value = keyValue[1].trim();

        value = removeArrayBrackets(value);
        String[] outputPairs = splitJsonElements(value);

        String irts=getSelectedElement(outputPairs, "{\"name\":\"ccs\",");
        irts = removeGroupingBrackets(irts);
        String[] irtPairs = splitJsonElements(irts);
        String irtdata=getSelectedElement(irtPairs, "\"data\":");
        String[] irtdataKeyValue = irtdata.split(":", 2);
        String irtdataValue = removeArrayBrackets(irtdataKeyValue[1]);
        TFloatArrayList irtdataArrayList=new TFloatArrayList();
        for (String num : irtdataValue.split(",")) {
        	irtdataArrayList.add(Float.parseFloat(num));
		}
        for (int i = 0; i < peptides.size(); i++) {
        	KoinaPrecursor precursor=peptides.get(i);
        	precursor.setIMS(irtdataArrayList.get(i));
		}
    }

    private static void parseRTData(String json, List<KoinaPrecursor> peptides) {
        json=removeGroupingBrackets(json);
        String[] pairs = splitJsonElements(json);
        
        String outputs=getSelectedElement(pairs, "\"outputs\":");
        String[] keyValue = outputs.split(":", 2);
        String value = keyValue[1].trim();

        value = removeArrayBrackets(value);
        String[] outputPairs = splitJsonElements(value);

        String irts=getSelectedElement(outputPairs, "{\"name\":\"irt\",");
        irts = removeGroupingBrackets(irts);
        String[] irtPairs = splitJsonElements(irts);
        String irtdata=getSelectedElement(irtPairs, "\"data\":");
        String[] irtdataKeyValue = irtdata.split(":", 2);
        String irtdataValue = removeArrayBrackets(irtdataKeyValue[1]);
        TFloatArrayList irtdataArrayList=new TFloatArrayList();
        for (String num : irtdataValue.split(",")) {
        	irtdataArrayList.add(Float.parseFloat(num));
		}
        for (int i = 0; i < peptides.size(); i++) {
        	KoinaPrecursor precursor=peptides.get(i);
        	precursor.setiRT(irtdataArrayList.get(i));
		}
    }

    private static void parseFragmentationData(String json, List<KoinaPrecursor> peptides) {
        json=removeGroupingBrackets(json);
        String[] pairs = splitJsonElements(json);
        
        String outputs=getSelectedElement(pairs, "\"outputs\":");
        String[] keyValue = outputs.split(":", 2);
        String value = keyValue[1].trim();

        value = removeArrayBrackets(value);
        String[] outputPairs = splitJsonElements(value);
        
        String mz=getSelectedElement(outputPairs, "{\"name\":\"mz\",");
        mz = removeGroupingBrackets(mz);
        String[] mzPairs = splitJsonElements(mz);

        String shape=getSelectedElement(mzPairs, "\"shape\":");
        String[] shapeKeyValue = shape.split(":", 2);
        String shapeValue = removeArrayBrackets(shapeKeyValue[1]);
        int n=Integer.parseInt(shapeValue.split(",")[0]);
        int m=Integer.parseInt(shapeValue.split(",")[1]);

        String mzdata=getSelectedElement(mzPairs, "\"data\":");
        String[] mzdataKeyValue = mzdata.split(":", 2);
        String mzdataValue = removeArrayBrackets(mzdataKeyValue[1]);
        TDoubleArrayList mzdataArrayList=new TDoubleArrayList();
        for (String num : mzdataValue.split(",")) {
			mzdataArrayList.add(Double.parseDouble(num));
		}

        String intensities=getSelectedElement(outputPairs, "{\"name\":\"intensities\",");
        intensities = removeGroupingBrackets(intensities);
        String[] intensitiesPairs = splitJsonElements(intensities);
        String intdata=getSelectedElement(intensitiesPairs, "\"data\":");
        String[] intdataKeyValue = intdata.split(":", 2);
        String intdataValue = removeArrayBrackets(intdataKeyValue[1]);
        TFloatArrayList intdataArrayList=new TFloatArrayList();
        for (String num : intdataValue.split(",")) {
			intdataArrayList.add(Float.parseFloat(num));
		}

        double[][] mzReshape=MatrixMath.reshape(mzdataArrayList.toArray(), n, m);
        float[][] intReshape=MatrixMath.reshape(intdataArrayList.toArray(), n, m);
        for (int i = 0; i < peptides.size(); i++) {
        	KoinaPrecursor precursor=peptides.get(i);
        	precursor.setIntensities(intReshape[i]);
        	precursor.setMzs(mzReshape[i]);
		}
    }

	private static String removeGroupingBrackets(String value) {
		value = value.trim();
		if (value.startsWith("{") && value.endsWith("}")) {
        	value = value.substring(1, value.length() - 1);
        }
		return value;
	}

	private static String removeArrayBrackets(String value) {
		value = value.trim();
		if (value.startsWith("[") && value.endsWith("]")) {
        	value = value.substring(1, value.length() - 1);
        }
		return value;
	}
    
    private static String getSelectedElement(String[] pairs, String key) {
        for (int i = 0; i < pairs.length; i++) {
        	if (pairs[i].startsWith(key)) {
        		return pairs[i];
        	}
        }
        return null;
    }

    private static String[] splitJsonElements(String json) {
        List<String> elements = new ArrayList<>();
        int braceLevel = 0;
        int bracketLevel = 0;
        int start = 0;

        for (int i = 0; i < json.length(); i++) {
            char ch = json.charAt(i);
            if (ch == '{') braceLevel++;
            else if (ch == '}') braceLevel--;
            else if (ch == '[') bracketLevel++;
            else if (ch == ']') bracketLevel--;
            else if (ch == ',' && braceLevel == 0 && bracketLevel == 0) {
                elements.add(json.substring(start, i));
                start = i + 1;
            }
        }
        elements.add(json.substring(start));

        return elements.toArray(new String[0]);
    }
	
	private static class KoinaPrecursor {
		private final String prositSequence;
		private final float nce;
		private final byte charge;
		private volatile float iRT=Float.NEGATIVE_INFINITY;
		private volatile float IMS=Float.NEGATIVE_INFINITY;
		private volatile float[] intensities=null;
		private volatile double[] mzs=null;
		HashSet<String> accessions=new HashSet<String>();
		
		public KoinaPrecursor(String peptideSequence, float nce, byte charge) {
			this.prositSequence=peptideSequence.replace("C", "C[UNIMOD:4]");
			this.nce = nce;
			this.charge = charge;
		}
		
		public void addAccession(String accession) {
			accessions.add(accession);
		}

		public void addAccessions(HashSet<String> multipleaccessions) {
			accessions.addAll(multipleaccessions);
		}
		
		public void setIMS(float IMS) {
			this.IMS = IMS;
		}
		
		public void setiRT(float iRT) {
			this.iRT = iRT;
		}
		
		public void setIntensities(float[] intensities) {
			this.intensities = intensities;
		}
		
		public void setMzs(double[] mzs) {
			this.mzs = mzs;
		}
		
		@Override
		public int hashCode() {
			return prositSequence.hashCode()+16807*Float.hashCode(nce)+1771561*Byte.hashCode(charge);
		}
		
		@Override
		public boolean equals(Object obj) {
			if (obj==null||!(obj instanceof KoinaPrecursor)) {
				return false;
			}
			KoinaPrecursor k=(KoinaPrecursor)obj;
			if (!prositSequence.equals(k.prositSequence)) {
				return false;
			}
			if (nce!=k.nce) {
				return false;
			}
			if (charge!=k.charge) {
				return false;
			}
			return true;
		}
		
		public AnnotatedLibraryEntry toEntry(AminoAcidConstants constants, SearchParameters params) {
			String peptideModSeq=prositSequence.replace("C[UNIMOD:4]", "C[+57.0]");
			double precursorMZ=constants.getChargedMass(peptideModSeq, charge);
			
			LibraryEntry e;
			if (Float.isInfinite(IMS)) {
				e=new LibraryEntry("Prosit", accessions, precursorMZ, charge, peptideModSeq, 1, iRT * 60f,
						0.0f, mzs, intensities, Optional.empty(), constants);
			} else {
				e=new LibraryEntry("Prosit", accessions, precursorMZ, charge, peptideModSeq, 1, iRT * 60f,
						0.0f, mzs, intensities, Optional.of(IMS), constants);
			}
			return new AnnotatedLibraryEntry(e, params);
		}
	}
}
