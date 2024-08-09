package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import org.deeplearning4j.datasets.iterator.utilty.ListDataSetIterator;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.layers.Convolution1DLayer;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.DropoutLayer;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
import org.deeplearning4j.nn.conf.layers.ZeroPadding1DLayer;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.nn.weights.WeightInit;
import org.deeplearning4j.optimize.listeners.ScoreIterationListener;
import org.deeplearning4j.util.ModelSerializer;
import org.nd4j.linalg.activations.Activation;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;
import org.nd4j.linalg.learning.config.Adam;
import org.nd4j.linalg.lossfunctions.LossFunctions;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class PeptidePropertiesPredictor {
	private static final int REPORTING_BATCH_SIZE = 1000;
	private static final int INITIAL_BATCH_SIZE = 64;
	private static final int MAX_BATCH_SIZE = 1024;
	private static final int EPOCHS_TO_2X_BATCH = 3;
    public static final int EMBED_DIMENSION = 128;
    public static final int N_RESNET_BLOCKS = 3;
    public static final int KERNEL_SIZE = 9;
    public static final int NUM_EPOCHS = 10;

    public static void main(String[] args) throws Exception {
    	File dir=new File("/Users/searleb/Documents/encyclopedia/prosit_examples_final/");
    	File modelLocation=new File(dir, "peptide_prediction_model.dl4j");
    	
    	SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
    	AminoAcidConstants aaConstants=parameters.getAAConstants();

		HashMap<String, DataSet> dataSetBySequence=new HashMap<String, DataSet>();
    	File[] listFiles = dir.listFiles(new SimpleFilenameFilter(".dlib"));
    	
    	// FIXME
    	//listFiles=new File[] {new File(dir, "UP000005640_9606.fasta.lys-n.z3_nce33.dlib.z3_nce33.dlib")};
    	
		for (File f : listFiles) {
			if (f.getName().indexOf(".trypsin")<0) {
				continue;
			}
    		String rtName=f.getName().substring(0, f.getName().length()-".z3_nce33.dlib".length())+".txt_rts.txt";
    		Logger.logLine("Reading RT file: "+rtName);

        	TObjectFloatHashMap<String> rtByPeptideModSeq=new TObjectFloatHashMap<String>();
    		TableParser.parseTSV(new File(f.getParent(), rtName), new TableParserMuscle() {
				
				@Override
				public void processRow(Map<String, String> row) {
					String peptideModSeq=row.get("PeptideModSeq");
					float rtInSec=Float.parseFloat(row.get("Pred_HI"))*60f;
					rtByPeptideModSeq.put(peptideModSeq, rtInSec);
				}
				
				@Override
				public void cleanup() {
				}
			});
    		
    		Logger.logLine("Reading Library file: "+f.getName());
    		LibraryFile library=new LibraryFile();
    		library.openFile(f);
    		ArrayList<LibraryEntry> entries=library.getAllEntries(false, aaConstants);
    		for (LibraryEntry entry : entries) {
				String key=entry.getPeptideModSeq()+"_"+entry.getPrecursorCharge();
				if (rtByPeptideModSeq.contains(entry.getPeptideModSeq())) {
					PeptideEncoding encoding=new PeptideEncoding(entry, rtByPeptideModSeq.get(entry.getPeptideModSeq()), parameters);
					dataSetBySequence.put(key, encoding.encodeDataset(parameters));
				} else {
					Logger.errorLine("Missing retention time for "+key);
				}
			}
    		library.close();
    	}

    	Logger.logLine("Loaded data "+dataSetBySequence.size()+" total precursors.");

        MultiLayerNetwork model=createModel(PeptideEncoding.ENCODED_INPUT_SIZE, PeptideEncoding.ENCODED_OUTPUT_SIZE, EMBED_DIMENSION, KERNEL_SIZE, N_RESNET_BLOCKS);
        Logger.logLine(model.summary());
        
        int batchSize=INITIAL_BATCH_SIZE;
    	double bestMeanAbsoluteError=Double.MAX_VALUE;
    	
    	ArrayList<DataSet> dataSet=new ArrayList<DataSet>(dataSetBySequence.values());
        for(int e=1; e<=NUM_EPOCHS; e++ ){
        	if (e%EPOCHS_TO_2X_BATCH==0) {
        		batchSize=batchSize*2;
        	}
        	Logger.logLine("Starting epoch "+e+", batch size: "+batchSize);
        	Collections.shuffle(dataSet);

            DataSetIterator iterator = new ListDataSetIterator<>(dataSet, batchSize);
            model.fit(iterator);
            
            int totalCount=0;
            double meanAbsoluteError=0.0;
            DataSetIterator fullIterator = new ListDataSetIterator<>(dataSet, MAX_BATCH_SIZE);
            INDArray results=model.output(fullIterator);

            TFloatArrayList deltas=new TFloatArrayList();
            for (int i = 0; i < dataSet.size(); i++) {
				double prediction=results.getDouble(i);
				double actual=dataSet.get(i).getLabels().getDouble(0);
				deltas.add((float)(actual-prediction));
				meanAbsoluteError+=Math.abs(actual-prediction);
			}
            totalCount+=deltas.size();
        	
            meanAbsoluteError=meanAbsoluteError/totalCount;
        	Logger.logLine("Epoch "+e+" mean absolute error: "+meanAbsoluteError+" from "+totalCount+" total peptides");

        	if (meanAbsoluteError<bestMeanAbsoluteError) {
        		bestMeanAbsoluteError=meanAbsoluteError;
	            try {
	            	ModelSerializer.writeModel(model, modelLocation, true);
	            	Logger.logLine("Model saved successfully.");
	            } catch (IOException ioe) {
	            	ioe.printStackTrace();
	            }
        	} else {
            	Logger.logLine("Final model wasn't better, so skipping temporary save.");
        	}
        }
    }
    
    private static MultiLayerNetwork createModel(int inputSize, int outputSize, int embedDim, int kernelSize, int resnetBlocks) {
    	NeuralNetConfiguration.ListBuilder builder = new NeuralNetConfiguration.Builder()
                .updater(new Adam(0.001))
                .dataType(PeptideEncoding.DEFAULT_DATA_TYPE)
                .list();

        // embedding layer
        builder.layer(new DenseLayer.Builder()
                .nIn(inputSize)
                .nOut(embedDim)
                .activation(Activation.RELU)
                .weightInit(WeightInit.XAVIER)
                .build());
        
        // resnet blocks
        for (int dilationRate = 1; dilationRate <=resnetBlocks; dilationRate++) {
            builder.layer(new Convolution1DLayer.Builder()
    		        .kernelSize(1)
    		        .stride(1)
    		        .nIn(embedDim)
    		        .nOut(embedDim)
    		        .activation(Activation.RELU)
    		        .weightInit(WeightInit.XAVIER)
    		        .build());
            
            builder.layer(new ZeroPadding1DLayer.Builder(dilationRate * (kernelSize - 1) / 2).build());

            builder.layer(new Convolution1DLayer.Builder()
    		        .kernelSize(kernelSize)
    		        .stride(1)
    		        .dilation(dilationRate)
    		        .nOut(embedDim)
    		        .activation(Activation.RELU)
    		        .weightInit(WeightInit.XAVIER)
    		        .build());
		}
        
        // dropout layer (keep 90% of changes)
        builder.layer(new DropoutLayer.Builder(0.9).build());
        
        // output layer
        builder.layer(new OutputLayer.Builder(LossFunctions.LossFunction.MSE)
                .nIn(embedDim)
                .nOut(outputSize)
                .activation(Activation.IDENTITY)
                .weightInit(WeightInit.XAVIER)
                .build());

        MultiLayerNetwork model = new MultiLayerNetwork(builder.build());
        model.init();
        model.setListeners(new ScoreIterationListener(REPORTING_BATCH_SIZE));

        return model;
    }
}