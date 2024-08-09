package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.util.ModelSerializer;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursorWithProteins;
import edu.washington.gs.maccoss.encyclopedia.gui.dia.WindowingSchemeWizard;

public class PeptidePropertiesModel {
	private static PeptidePropertiesModel staticModel=null; 
	private final MultiLayerNetwork model;
	
	public static synchronized PeptidePropertiesModel getModel() throws IOException {
		if (staticModel==null) {
			InputStream is=WindowingSchemeWizard.class.getResourceAsStream("/models/peptide_prediction_model.dl4j");
	        MultiLayerNetwork model = ModelSerializer.restoreMultiLayerNetwork(is);
	        
	        staticModel=new PeptidePropertiesModel(model);
		}
		return staticModel;
	}

	private PeptidePropertiesModel(MultiLayerNetwork model) {
		this.model = model;
	}
	
	public ArrayList<LibraryEntry> predict(ArrayList<PeptidePrecursorWithProteins> peptides, AminoAcidConstants constants) {
		ArrayList<PeptidePrecursorWithProteins> predictedPrecursors=new ArrayList<PeptidePrecursorWithProteins>();
        ArrayList<INDArray> encodings=new ArrayList<INDArray>();
        for (PeptidePrecursorWithProteins peptide : peptides) {
        	if (peptide.getPeptideSeq().length()>PeptideEncoding.MAX_PEPTIDE_LENGTH-2||peptide.getPrecursorCharge()>PeptideEncoding.MAX_CHARGE) {
        		continue;
        	}
        	predictedPrecursors.add(peptide);
			encodings.add(PeptideEncoding.encodeInput(peptide.getPeptideModSeq(), peptide.getPrecursorCharge(), constants));
		}
        
        INDArray inputData = Nd4j.concat(0, encodings.toArray(new INDArray[0]));
        inputData = inputData.reshape(encodings.size(), PeptideEncoding.ENCODED_INPUT_SIZE);
        INDArray outputData = model.output(inputData);
        outputData=outputData.reshape(encodings.size(), PeptideEncoding.ENCODED_OUTPUT_SIZE);
        
        ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
        for (int i = 0; i < outputData.rows(); i++) {
			INDArray row=outputData.getRow(i);
			PeptidePrecursorWithProteins peptide=predictedPrecursors.get(i);
			entries.add(PeptideEncoding.outputToEntry(peptide.getPeptideModSeq(), peptide.getPrecursorCharge(), peptide.getAccessions(), row, constants));
		}
        return entries;
	}

}
