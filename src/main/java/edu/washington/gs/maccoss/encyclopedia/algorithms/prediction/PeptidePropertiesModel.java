package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.util.ModelSerializer;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.DataSetUtil;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;
import org.nd4j.linalg.factory.Nd4j;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursorWithProteins;
import edu.washington.gs.maccoss.encyclopedia.gui.dia.WindowingSchemeWizard;

public class PeptidePropertiesModel {
	private static PeptidePropertiesModel staticModel=null; 
	private final ComputationGraph model;
	
	public static synchronized PeptidePropertiesModel getModel() throws IOException {
		if (staticModel==null) {
			InputStream is=WindowingSchemeWizard.class.getResourceAsStream("/models/peptide_prediction_model.dl4j");
			ComputationGraph model = ModelSerializer.restoreComputationGraph(is);
	        
	        staticModel=new PeptidePropertiesModel(model);
		}
		return staticModel;
	}

	private PeptidePropertiesModel(ComputationGraph model) {
		this.model = model;
	}
	
	public ArrayList<LibraryEntry> predict(ArrayList<PeptidePrecursorWithProteins> peptides, AminoAcidConstants constants) {
		ArrayList<LibraryEntry> predictedPrecursors=new ArrayList<LibraryEntry>();
		
		int batchSize=PeptidePropertiesPredictor.MAX_BATCH_SIZE;
		
		ArrayList<PeptidePrecursorWithProteins> encodedPeptides=new ArrayList<PeptidePrecursorWithProteins>();
		ArrayList<INDArray[]> encodings=new ArrayList<INDArray[]>();
        for (PeptidePrecursorWithProteins peptide : peptides) {
        	if (peptide.getPeptideSeq().length()>PeptideEncoding.MAX_PEPTIDE_LENGTH-2||peptide.getPrecursorCharge()>PeptideEncoding.MAX_CHARGE) {
        		continue;
        	}
        	encodedPeptides.add(peptide);
			encodings.add(PeptideEncoding.encodeInput(peptide.getPeptideModSeq(), peptide.getPrecursorCharge(), constants));
			
			if (encodings.size()>=batchSize) {
				// run batch
		        INDArray[][] arr = encodings.toArray(new INDArray[encodings.size()][0]);
		        INDArray[] input=DataSetUtil.mergeFeatures(arr, null).getFirst();
		        INDArray[] output=model.output(input);
		        
		        INDArray fragmentationResults=output[0];
		        INDArray rtResults=output[1];
		        INDArray imsResults=output[2];
		        

		        for (int i = 0; i < encodings.size(); i++) {
		        	float[] fragmentData=new float[PeptideEncoding.ENCODED_OUTPUT_FRAGMENT_SIZE];
		        	int start=i*PeptideEncoding.ENCODED_OUTPUT_FRAGMENT_SIZE;
		        	for (int j = 0; j < fragmentData.length; j++) {
		        		fragmentData[j]=fragmentationResults.getFloat(start+j);
					}
		        	
		        	float rt=rtResults.getFloat(i);
		        	float ims=imsResults.getFloat(i);
		        	PeptidePrecursorWithProteins thisPeptide=encodedPeptides.get(i);
		        	LibraryEntry entry=PeptideEncoding.outputToEntry(thisPeptide.getPeptideModSeq(), fragmentData, rt, ims, thisPeptide.getPrecursorCharge(), thisPeptide.getAccessions(), constants);
		        	
		        	predictedPrecursors.add(entry);
		        }
				
				encodings.clear();
				encodedPeptides.clear();
			}
		}
        return predictedPrecursors;
	}

}
