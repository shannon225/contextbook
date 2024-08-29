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

import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.MatrixMath;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;

public abstract class PrositFragmentationPredictionModel implements KoinaFeaturePredictionModel {

    public abstract boolean useNCE();

	public void updateResults(List<KoinaPrecursor> peptides, ArrayList<Pair<double[], float[]>> values) {
		assert(peptides.size()==values.size());
		for (int i = 0; i < peptides.size(); i++) {
			Pair<double[], float[]> pair = values.get(i);
			peptides.get(i).setMzs(pair.getX());
			peptides.get(i).setIntensities(pair.getY());
		}
	}
	
	@Override
	public boolean canModelPeptide(AminoAcidEncoding[] aas, byte precursorCharge) {
		return CommonModelConstraints.canModelPeptidePrositStandard(aas, precursorCharge);
	}
    
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
	
	        JSONObject precursorCharges = new JSONObject();
	        precursorCharges.put("name", "precursor_charges");
	        precursorCharges.put("shape", new JSONArray("["+peptides.size()+",1]"));
	        precursorCharges.put("datatype", "INT32");
	        precursorCharges.put("data", new JSONArray(charges.toArray()));
	
	        inputsArray.put(peptideSequences);
	        
	        if (useNCE()) {
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

	        ArrayList<Pair<double[], float[]>> results=parseFragmentationData(response.toString());
	        updateResults(peptides, results);
	        
		} catch (IOException ioe) {
			throw new EncyclopediaException("IO error getting Koina result", ioe);
		}
	}
    
    public ArrayList<Pair<double[], float[]>> parseFragmentationData(String json) {
        json=JSONParsingUtilities.removeGroupingBrackets(json);
        String[] pairs = JSONParsingUtilities.splitJsonElements(json);
        
        String outputs=JSONParsingUtilities.getSelectedElement(pairs, "\"outputs\":");
        String[] keyValue = outputs.split(":", 2);
        String value = keyValue[1].trim();

        value = JSONParsingUtilities.removeArrayBrackets(value);
        String[] outputPairs = JSONParsingUtilities.splitJsonElements(value);
        
        String mz=JSONParsingUtilities.getSelectedElement(outputPairs, "{\"name\":\"mz\",");
        mz = JSONParsingUtilities.removeGroupingBrackets(mz);
        String[] mzPairs = JSONParsingUtilities.splitJsonElements(mz);

        String shape=JSONParsingUtilities.getSelectedElement(mzPairs, "\"shape\":");
        String[] shapeKeyValue = shape.split(":", 2);
        String shapeValue = JSONParsingUtilities.removeArrayBrackets(shapeKeyValue[1]);
        int n=Integer.parseInt(shapeValue.split(",")[0]);
        int m=Integer.parseInt(shapeValue.split(",")[1]);

        String mzdata=JSONParsingUtilities.getSelectedElement(mzPairs, "\"data\":");
        String[] mzdataKeyValue = mzdata.split(":", 2);
        String mzdataValue = JSONParsingUtilities.removeArrayBrackets(mzdataKeyValue[1]);
        TDoubleArrayList mzdataArrayList=new TDoubleArrayList();
        for (String num : mzdataValue.split(",")) {
			mzdataArrayList.add(Double.parseDouble(num));
		}

        String intensities=JSONParsingUtilities.getSelectedElement(outputPairs, "{\"name\":\"intensities\",");
        intensities = JSONParsingUtilities.removeGroupingBrackets(intensities);
        String[] intensitiesPairs = JSONParsingUtilities.splitJsonElements(intensities);
        String intdata=JSONParsingUtilities.getSelectedElement(intensitiesPairs, "\"data\":");
        String[] intdataKeyValue = intdata.split(":", 2);
        String intdataValue = JSONParsingUtilities.removeArrayBrackets(intdataKeyValue[1]);
        TFloatArrayList intdataArrayList=new TFloatArrayList();
        for (String num : intdataValue.split(",")) {
			intdataArrayList.add(Float.parseFloat(num));
		}

        double[][] mzReshape=MatrixMath.reshape(mzdataArrayList.toArray(), n, m);
        float[][] intReshape=MatrixMath.reshape(intdataArrayList.toArray(), n, m);
        ArrayList<Pair<double[], float[]>> patterns=new ArrayList<Pair<double[],float[]>>();
        for (int i = 0; i < mzReshape.length; i++) {
        	patterns.add(new Pair<double[], float[]>(mzReshape[i], intReshape[i]));
		}
        return patterns;
    }
}
