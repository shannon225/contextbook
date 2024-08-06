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

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class PeptidePropertiesPredictor {
	private static final int INITIAL_BATCH_SIZE = 64;
	private static final int MAX_BATCH_SIZE = 1024;
	private static final int EPOCHS_TO_2X_BATCH = 2;
    public static final int EMBED_DIMENSION = 128;
    public static final int N_RESNET_BLOCKS = 3;
    public static final int KERNEL_SIZE = 7;
    public static final int NUM_EPOCHS = 8;

    private static final int MAX_CHARGE=6;
    private static final int MAX_PEPTIDE_LENGTH=32; // number of termini + 30
    
    private static final int ENCODED_INPUT_SIZE=MAX_PEPTIDE_LENGTH*EncodedAminoAcid.MAX_ENCODING_LENGTH+MAX_CHARGE;
    private static final int ENCODED_OUTPUT_SIZE=MAX_PEPTIDE_LENGTH*4+2;
    
    private static final AminoAcidConstants aaConstants=new AminoAcidConstants();

    public static void main(String[] args) throws Exception {
    	File dir=new File("/Users/searleb/Documents/encyclopedia/prosit_examples_final/");
    	File modelLocation=new File(dir, "peptide_prediction_model.dl4j");
    	
    	SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
    	AminoAcidConstants aaConstants=parameters.getAAConstants();
    	
    	TObjectFloatHashMap<String> rtByPeptideModSeq=new TObjectFloatHashMap<String>();
    	for (File f : dir.listFiles(new SimpleFilenameFilter(".txt_rts.txt"))) {
    		Logger.logLine("Reading RT file: "+f.getName());
    		TableParser.parseTSV(f, new TableParserMuscle() {
				
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
    	}
    	
    	Logger.logLine("Read RTs for "+rtByPeptideModSeq.size()+" total peptides.");

		HashMap<String, DataSet> dataSetBySequence=new HashMap<String, DataSet>();
    	for (File f : dir.listFiles(new SimpleFilenameFilter(".dlib"))) {
    		Logger.logLine("Reading Library file: "+f.getName());
    		
    		LibraryFile library=new LibraryFile();
    		library.openFile(f);
    		ArrayList<LibraryEntry> entries=library.getAllEntries(false, aaConstants);
    		for (LibraryEntry entry : entries) {
				String key=entry.getPeptideModSeq()+"_"+entry.getPrecursorCharge();
				if (rtByPeptideModSeq.contains(entry.getPeptideModSeq())) {
					PeptideEncoding encoding=new PeptideEncoding(entry, rtByPeptideModSeq.get(entry.getPeptideModSeq()), parameters);
					dataSetBySequence.put(key, new DataSet(encoding.encodeInput(parameters), encoding.encodeResult()));
				} else {
					Logger.errorLine("Missing retention time for "+key);
				}
			}
    		library.close();
    	}

    	Logger.logLine("Loaded data "+dataSetBySequence.size()+" total precursors.");

        MultiLayerNetwork model=createModel(ENCODED_INPUT_SIZE, ENCODED_OUTPUT_SIZE, EMBED_DIMENSION, KERNEL_SIZE, N_RESNET_BLOCKS);
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
    
    public static void mainPartialModelTraining(String[] args) throws Exception {
    	File f=new File("/Users/searleb/Documents/encyclopedia/prosit_examples_final/UP000005640_9606.fasta.trypsin.z6_nce33.dlib.z3_nce33.dlib");
    	File saveLocation=new File(f.getParentFile(), "peptide_prediction_model.dl4j");
    	SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
    	LibraryFile library=new LibraryFile();
		library.openFile(f);
		ArrayList<LibraryEntry> entries=library.getAllEntries(false, aaConstants);
		
		HashMap<String, DataSet> dataSets=new HashMap<String, DataSet>();
		for (LibraryEntry entry : entries) {
			String key=entry.getPeptideModSeq()+"_"+entry.getPrecursorCharge();
			PeptideEncoding encoding=new PeptideEncoding(entry, entry.getScanStartTime(), parameters);

			dataSets.put(key, new DataSet(encoding.encodeInput(parameters), encoding.encodeResult()));
		}
		
    	Logger.logLine("Loaded "+dataSets.size()+" total peptides.");

        MultiLayerNetwork model=createModel(ENCODED_INPUT_SIZE, ENCODED_OUTPUT_SIZE, EMBED_DIMENSION, KERNEL_SIZE, N_RESNET_BLOCKS);
        Logger.logLine(model.summary());
        

    	int batchSize=INITIAL_BATCH_SIZE;
    	double bestMeanAbsoluteError=Double.MAX_VALUE;
    	
    	ArrayList<DataSet> dataSet=new ArrayList<DataSet>(dataSets.values());
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
    
    private static class PeptideEncoding {
    	private final String peptideModSeq;
    	private final byte charge;
    	private final float rtInSec;
    	private final float ims;
    	private final float[] bp1=new float[MAX_PEPTIDE_LENGTH];
    	private final float[] bp2=new float[MAX_PEPTIDE_LENGTH];
    	private final float[] yp1=new float[MAX_PEPTIDE_LENGTH];
    	private final float[] yp2=new float[MAX_PEPTIDE_LENGTH];
    	

    	public INDArray encodeInput(SearchParameters parameters) {
        	EncodedAminoAcid[] aas=EncodedAminoAcid.getAAs(peptideModSeq, parameters.getAAConstants());
        	
            INDArray encoded = Nd4j.zeros(MAX_PEPTIDE_LENGTH, EncodedAminoAcid.MAX_ENCODING_LENGTH);
            
            int start=aas[0].isNTerm()?0:1;
            for (int i = start; i < aas.length; i++) {
                encoded.putScalar(new int[]{i, aas[i].getIndex()}, 1.0);
            }
            
            INDArray encodedCharge=Nd4j.zeros(1, MAX_CHARGE);
            encodedCharge.putScalar(new int[] {0, charge-1}, 1.0);
            
            INDArray reshape = encoded.reshape(1, MAX_PEPTIDE_LENGTH * EncodedAminoAcid.MAX_ENCODING_LENGTH);
            
			return Nd4j.concat(1, reshape, encodedCharge);
    	}
    	
    	public INDArray encodeResult() {
    		INDArray bp1array=Nd4j.create(bp1);
    		INDArray bp2array=Nd4j.create(bp2);
    		INDArray yp1array=Nd4j.create(yp1);
    		INDArray yp2array=Nd4j.create(yp2);
    		INDArray scalar=Nd4j.create(new float[] {rtInSec, ims});
			
			INDArray encoded=Nd4j.concat(0, bp1array, bp2array, yp1array, yp2array, scalar);
			return encoded.reshape(1, encoded.length());
    	}
    	
    	public PeptideEncoding(LibraryEntry entry, float rtInSec, SearchParameters parameters) {
    		this.peptideModSeq=entry.getPeptideModSeq();
    		this.charge=entry.getPrecursorCharge();
    		this.rtInSec=rtInSec;
    		this.ims=entry.getIonMobility().get();
    		
			double[] massArray = entry.getMassArray();
			float[] intensityArray=entry.getIntensityArray();
    		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
    		for (FragmentIon fragmentIon : model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), true)) {
				byte fragcharge=IonType.getCharge(fragmentIon.getType());
				if (fragcharge>3||fragcharge<1) {
					continue;
				}
				IonType type=IonType.getCanonicalIonType(fragmentIon.getType());
				float[] array=null;
				if (type==IonType.b) {
					switch (fragcharge) {
						case 1: array=bp1; break;
						case 2: array=bp2; break;
						default: break;
					}
				} else if (type==IonType.y) {
					switch (fragcharge) {
						case 1: array=yp1; break;
						case 2: array=yp2; break;
						default: break;
					}
				}
				if (array==null) {
					continue;
				}
				
				float intensity=parameters.getFragmentTolerance().getMaxIntensity(massArray, intensityArray, fragmentIon.getMass());
				int fragindex=fragmentIon.getIndex();
				
				array[fragindex]=intensity;
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
        model.setListeners(new ScoreIterationListener(100));

        return model;
    }
}