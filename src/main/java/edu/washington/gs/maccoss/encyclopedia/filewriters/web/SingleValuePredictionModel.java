package edu.washington.gs.maccoss.encyclopedia.filewriters.web;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import org.json.JSONArray;
import org.json.JSONObject;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;

public abstract class SingleValuePredictionModel implements KoinaFeaturePredictionModel {

    public abstract void updateResults(List<KoinaPrecursor> peptides, float[] values);
    public abstract boolean useCharge();
    public abstract boolean useNCE();
    public abstract String getDataTypeName();
    
    @Override
    public String toString() {
    	return getName();
    }
    
    @Override
	public void updatePeptides(List<KoinaPrecursor> peptides) {
		ArrayList<String> pepseqs=new ArrayList<String>();
		TFloatArrayList NCEs=new TFloatArrayList();
		TIntArrayList charges=new TIntArrayList();
		for (KoinaPrecursor pep : peptides) {
			pepseqs.add(pep.getKoinaSequence());
			NCEs.add(pep.getNCE());
			charges.add(pep.getCharge());
		}
		
		try {
			URL url=getURL();
			
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
	        
	        if (useCharge()) {
		        JSONObject precursorCharges = new JSONObject();
		        precursorCharges.put("name", "precursor_charges");
		        precursorCharges.put("shape", new JSONArray("["+peptides.size()+",1]"));
		        precursorCharges.put("datatype", "INT32");
		        precursorCharges.put("data", new JSONArray(charges.toArray()));
		        inputsArray.put(precursorCharges);
	        }

	        if (useNCE()) {
		        JSONObject collisionEnergies = new JSONObject();
		        collisionEnergies.put("name", "collision_energies");
		        collisionEnergies.put("shape", new JSONArray("["+peptides.size()+",1]"));
		        collisionEnergies.put("datatype", "FP32");
		        collisionEnergies.put("data", new JSONArray(NCEs.toArray()));
	        	inputsArray.put(collisionEnergies);
	        }
	        
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

	        float[] results=parseData(response.toString());
	        updateResults(peptides, results);
	        
		} catch (IOException ioe) {
			throw new EncyclopediaException("IO error getting Koina result", ioe);
		}
	}

    public float[] parseData(String json) {
        json=JSONParsingUtilities.removeGroupingBrackets(json);
        String[] pairs = JSONParsingUtilities.splitJsonElements(json);
        
        String outputs=JSONParsingUtilities.getSelectedElement(pairs, "\"outputs\":");
        String[] keyValue = outputs.split(":", 2);
        String value = keyValue[1].trim();

        value = JSONParsingUtilities.removeArrayBrackets(value);
        String[] outputPairs = JSONParsingUtilities.splitJsonElements(value);

        String irts=JSONParsingUtilities.getSelectedElement(outputPairs, "{\"name\":\""+getDataTypeName()+"\",");
        irts = JSONParsingUtilities.removeGroupingBrackets(irts);
        String[] irtPairs = JSONParsingUtilities.splitJsonElements(irts);
        String irtdata=JSONParsingUtilities.getSelectedElement(irtPairs, "\"data\":");
        String[] irtdataKeyValue = irtdata.split(":", 2);
        String irtdataValue = JSONParsingUtilities.removeArrayBrackets(irtdataKeyValue[1]);
        TFloatArrayList irtdataArrayList=new TFloatArrayList();
        for (String num : irtdataValue.split(",")) {
        	irtdataArrayList.add(Float.parseFloat(num));
		}
        return irtdataArrayList.toArray();
    }
}
