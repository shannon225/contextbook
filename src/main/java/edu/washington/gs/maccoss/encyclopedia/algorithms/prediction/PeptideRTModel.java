package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.util.ModelSerializer;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.gui.dia.WindowingSchemeWizard;

public class PeptideRTModel {
	private static PeptideRTModel staticModel=null; 
	private final MultiLayerNetwork model;
	
	public static synchronized PeptideRTModel getModel() throws IOException {
		if (staticModel==null) {
			InputStream is=WindowingSchemeWizard.class.getResourceAsStream("/models/peptide_rt_model.dl4j");
            MultiLayerNetwork model = ModelSerializer.restoreMultiLayerNetwork(is);
	        
	        staticModel=new PeptideRTModel(model);
		}
		return staticModel;
	}

	private PeptideRTModel(MultiLayerNetwork model) {
		this.model = model;
	}
	
	public ArrayList<LibraryEntry> predict(ArrayList<LibraryEntry> peptides, AminoAcidConstants constants) {
		ArrayList<LibraryEntry> predictedPrecursors=new ArrayList<>();
    	ArrayList<PeptideRTEncoding> encodings=new ArrayList<>();
		for (LibraryEntry entry : peptides) {
			AminoAcidEncoding[] aas=AminoAcidEncoding.getAAs(entry.getPeptideModSeq(), constants);
	    	if (aas.length>PeptideRTEncoding.MAXIMUM_PEPTIDE_LENGTH) {
	    		continue;
	    	}
	    	
	    	predictedPrecursors.add(entry);
			encodings.add(new PeptideRTEncoding(aas, entry.getRetentionTime()));
		}

        DataSetIterator iterator = new EncodedDataSetIterator(encodings, PeptideRTPredictor.MAX_BATCH_SIZE);
        INDArray output = model.output(iterator);
        
        for (int i = 0; i < predictedPrecursors.size(); i++) {
        	double prediction = output.getDouble(i);
        	predictedPrecursors.set(i, predictedPrecursors.get(i).updateRetentionTime((float)prediction));
        }
        return predictedPrecursors;
	}
}
