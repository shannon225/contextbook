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
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.learning.config.Adam;
import org.nd4j.linalg.lossfunctions.LossFunctions;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.TwoDimensionalKDE;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import gnu.trove.list.array.TFloatArrayList;

public class PeptideRTPredictor {
	private static final int INITIAL_BATCH_SIZE = 64;
	private static final int MAX_BATCH_SIZE = 1024;
	private static final int EPOCHS_TO_2X_BATCH = 3;
    public static final int EMBED_DIMENSION = 64;
    public static final int N_RESNET_BLOCKS = 3;
    public static final int KERNEL_SIZE = 7;
    public static final int NUM_EPOCHS = 10;
    
    private static final AminoAcidConstants aaConstants=new AminoAcidConstants();

    public static void main(String[] args) {
    	File db=new File("C:\\Users\\searl\\Documents\\projects\\Chronologer_DB_220308.txt");
    	db=new File("C:\\Users\\searl\\Downloads\\combined.txt\\combined.txt");
    	File saveLocation=new File(db.getParentFile(), "peptide_rt_model.dl4j");
    	
    	ArrayList<XYPoint> points=new ArrayList<XYPoint>();
    	
    	double meanAbsoluteError=0.0;
    	try {
    		LibraryFile library=new LibraryFile();
    		library.openFile(new File(db.getParentFile(), "23aug2017_hela_serum_timecourse_pool_wide_001_170829031834.dia.encyclopedia2.txt.elib"));
    		ArrayList<LibraryEntry> entries=library.getAllEntries(false, aaConstants);

            MultiLayerNetwork model = ModelSerializer.restoreMultiLayerNetwork(saveLocation);
            System.out.println("Model loaded successfully.");
            
            long time=System.currentTimeMillis();
            
            ArrayList<INDArray> encodings=new ArrayList<INDArray>();
            
            for (LibraryEntry entry : entries) {
				encodings.add(AminoAcidEncoding.encode(entry.getPeptideModSeq(), aaConstants, PeptideRTEncoding.MAXIMUM_PEPTIDE_LENGTH));
			}
            
            INDArray inputData = Nd4j.concat(0, encodings.toArray(new INDArray[0]));
            inputData = inputData.reshape(entries.size(), PeptideRTEncoding.MAXIMUM_PEPTIDE_LENGTH * AminoAcidEncoding.MAX_ENCODING_LENGTH);
            INDArray output = model.output(inputData);

            for (int i = 0; i < entries.size(); i++) {
            	points.add(new XYPoint(output.getDouble(i), entries.get(i).getScanStartTime()));
                //System.out.println("Prediction for sequence " + entries.get(i).getPeptideModSeq() + "\t" + output.getDouble(i)+ "\t" + entries.get(i).getScanStartTime());
            }
            
    		TwoDimensionalKDE twoDimKDE=new TwoDimensionalKDE(points, 3000);
    		Function rtWarper=twoDimKDE.trace();
            for (int i = 0; i < entries.size(); i++) {
            	meanAbsoluteError+=Math.abs(output.getDouble(i)-rtWarper.getXValue(entries.get(i).getScanStartTime()));
            	//System.out.println(entries.get(i).getPeptideModSeq()+"\t"+(float)output.getDouble(i)+"\t"+rtWarper.getXValue(entries.get(i).getScanStartTime()));
            }
            meanAbsoluteError=meanAbsoluteError/entries.size();

            System.out.println("time per prediction: "+(System.currentTimeMillis()-time)+" msec for "+entries.size()+" peptides");
            System.out.println("testing reserved "+meanAbsoluteError+" mean absolute error");
            
    	} catch (Exception e) {
    		e.printStackTrace();
    	}
    }

    public static void main2(String[] args) {
    	File db=new File("C:\\Users\\searl\\Documents\\projects\\Chronologer_DB_220308.txt");
    	db=new File("C:\\Users\\searl\\Downloads\\combined.txt\\combined.txt");
    	File saveLocation=new File(db.getParentFile(), "peptide_rt_model.dl4j");
    	
    	ArrayList<ArrayList<PeptideRTEncoding>> dataSets = new ArrayList<>();
    	ArrayList<String> dataSetNames = new ArrayList<>();
    	
    	TableParser.parseTSV(db, new TableParserMuscle() {
        	ArrayList<PeptideRTEncoding> currentDataSet=null;
        	String currentDataSetName=null;
        	int count=0;
        	
			@Override
			public void processRow(Map<String, String> row) {
				String source=row.get("Source");
				if (!source.equals(currentDataSetName)) {
					currentDataSet=new ArrayList<PeptideRTEncoding>();
					dataSets.add(currentDataSet);
					dataSetNames.add(source);
					currentDataSetName=source;
				}
				
				String peptideModSeq=row.get("PeptideModSeq");
				float hirt=Float.parseFloat(row.get("HI"));
				
				AminoAcidEncoding[] aas=AminoAcidEncoding.getAAs(peptideModSeq, aaConstants);
		    	if (aas.length>PeptideRTEncoding.MAXIMUM_PEPTIDE_LENGTH) return;
		    	
				currentDataSet.add(new PeptideRTEncoding(aas, hirt));

				count++;
				if (count%1000000==0) {
					Logger.logLine("Loading "+count+"...");
				}
			}
			
			@Override
			public void cleanup() {
			}
		});
    	
    	int count=0;
    	for (ArrayList<PeptideRTEncoding> dataset : dataSets) {
			count+=dataset.size();
		}
    	Logger.logLine("Loaded "+dataSets.size()+" datasets with "+count+" total peptides.");

        MultiLayerNetwork model=createModel(PeptideRTEncoding.MAXIMUM_PEPTIDE_LENGTH * AminoAcidEncoding.MAX_ENCODING_LENGTH, 1, 
        		EMBED_DIMENSION, KERNEL_SIZE, N_RESNET_BLOCKS);
        Logger.logLine(model.summary());

		ArrayList<PeptideRTEncoding> trimmedDataSet = new ArrayList<PeptideRTEncoding>();
    	
    	// start with full dataset
    	for (int i = 0; i < dataSets.size(); i++) {
    		trimmedDataSet.addAll(dataSets.get(i));
		}
    	
    	int batchSize=INITIAL_BATCH_SIZE;
    	double bestMeanAbsoluteError=Double.MAX_VALUE;
    	
        for(int e=1; e<=NUM_EPOCHS; e++ ){
        	if (e%EPOCHS_TO_2X_BATCH==0) {
        		batchSize=batchSize*2;
        	}
        	Logger.logLine("Starting epoch "+e+", batch size: "+batchSize);
        	Collections.shuffle(trimmedDataSet);

            DataSetIterator iterator = new EncodedDataSetIterator(trimmedDataSet, batchSize);
            model.fit(iterator);

            trimmedDataSet.clear();

            int totalCount=0;
            double meanAbsoluteError=0.0;
        	for (int d = 0; d < dataSets.size(); d++) {
        		ArrayList<PeptideRTEncoding> dataSet=dataSets.get(d);
	            DataSetIterator fullIterator = new EncodedDataSetIterator(dataSet, MAX_BATCH_SIZE);
	            INDArray results=model.output(fullIterator);

	            TFloatArrayList deltas=new TFloatArrayList();
                for (int i = 0; i < dataSet.size(); i++) {
					double prediction=results.getDouble(i);
					float actual=dataSet.get(i).getHiRT();
					deltas.add(actual-(float)prediction);
					meanAbsoluteError+=Math.abs(actual-prediction);
				}
                totalCount+=deltas.size();
                
                float[] deltaArray=deltas.toArray();
                float low=QuickMedian.select(deltaArray, 0.005f);
                float high=QuickMedian.select(deltaArray, 0.995f);
                for (int i = 0; i < dataSet.size(); i++) {
					double prediction=results.getDouble(i);
					float actual=dataSet.get(i).getHiRT();
					double delta=actual-prediction;
					if (delta<high&&delta>low) {
						trimmedDataSet.add(dataSet.get(i));
					}
				}
        	}
        	
            meanAbsoluteError=meanAbsoluteError/totalCount;
        	Logger.logLine("Epoch "+e+" mean absolute error: "+meanAbsoluteError+" from "+totalCount+" total peptides");

        	if (meanAbsoluteError<bestMeanAbsoluteError) {
        		bestMeanAbsoluteError=meanAbsoluteError;
	            try {
	            	ModelSerializer.writeModel(model, saveLocation, true);
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
        model.setListeners(new ScoreIterationListener(1000));

        return model;
    }
}