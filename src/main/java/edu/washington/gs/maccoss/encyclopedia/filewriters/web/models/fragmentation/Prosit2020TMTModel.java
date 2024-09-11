package edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.json.JSONArray;
import org.json.JSONObject;

import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.CommonModelConstraints;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaPrecursor;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.PrositFragmentationPredictionModel;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;

public class Prosit2020TMTModel extends PrositFragmentationPredictionModel {
	private static final String CID="CID";
	private static final String HCD="HCD";
	
	private final boolean isHCD;
	
	public Prosit2020TMTModel(boolean isHCD) {
		super();
		this.isHCD = isHCD;
	}
	
	@Override
	public boolean canModelPeptide(AminoAcidEncoding[] aas, byte precursorCharge) {
		return CommonModelConstraints.canModelPeptidePrositTMT(aas, precursorCharge);
	}

	@Override
	public String getName() {
		return "Prosit 2020 TMT "+(isHCD?HCD:CID);
	}

	@Override
	public URL getURL(String baseURL) {
		try {
			return new URL(baseURL+"v2/models/Prosit_2020_intensity_TMT/infer");
		} catch (MalformedURLException e) {
			throw new EncyclopediaException("Error getting Koina URL", e);
		}
	}

	@Override
	public boolean useNCE() {
		return false;
	}

    @Override
    public void updatePeptides(List<KoinaPrecursor> peptides, String baseURL) {
		ArrayList<String> pepseqs=new ArrayList<String>();
		TFloatArrayList NCEs=new TFloatArrayList();
		TIntArrayList charges=new TIntArrayList();
		for (KoinaPrecursor pep : peptides) {
			pepseqs.add(pep.getKoinaSequence());
			NCEs.add(pep.getNCE());
			charges.add(pep.getCharge());
		}
		
		try {
			URL url=getURL(baseURL);
			
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
	        inputsArray.put(precursorCharges);
	        
	        // must use NCE
	        JSONObject collisionEnergies = new JSONObject();
	        collisionEnergies.put("name", "collision_energies");
	        collisionEnergies.put("shape", new JSONArray("["+peptides.size()+",1]"));
	        collisionEnergies.put("datatype", "FP32");
	        collisionEnergies.put("data", new JSONArray(NCEs.toArray()));
        	inputsArray.put(collisionEnergies);

	        String[] fragTypeArray=new String[peptides.size()];
	        if (isHCD) {
	        	Arrays.fill(fragTypeArray, HCD);
	        } else {
	        	Arrays.fill(fragTypeArray, CID);
	        }
	        JSONObject fragType = new JSONObject();
	        fragType.put("name", "fragmentation_types");
	        fragType.put("shape", new JSONArray("["+peptides.size()+",1]"));
	        fragType.put("datatype", "BYTES");
	        fragType.put("data", new JSONArray(fragTypeArray));
	        inputsArray.put(fragType);
	
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

	        ArrayList<Pair<double[], float[]>> results=parseFragmentationData(response.toString());
	        updateResults(peptides, results);
	        
		} catch (IOException ioe) {
			throw new EncyclopediaException("IO error getting Koina result", ioe);
		}
	}

}
