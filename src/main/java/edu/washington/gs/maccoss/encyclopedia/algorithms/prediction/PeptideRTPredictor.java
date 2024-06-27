package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
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
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.PivotTableGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import gnu.trove.list.array.TFloatArrayList;

public class PeptideRTPredictor {
	private static final int MAX_BATCH_SIZE = 1024;
	private static final String AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY";
    private static final int AMINO_ACID_COUNT = AMINO_ACIDS.length();
    private static final int MAXIMUM_PEPTIDE_LENGTH = 35;
    public static final int EMBED_DIMENSION = 64;
    public static final int N_RESNET_BLOCKS = 3;
    public static final int KERNEL_SIZE = 7;
    public static final int NUM_EPOCHS = 30;
    
    private static final AminoAcidConstants aaConstants=new AminoAcidConstants();

    public static void main(String[] args) {
    	File db=new File("/Users/searleb/Downloads/Chronologer_DB_220308.txt");
    	File saveLocation=new File(db.getParentFile(), "peptide_rt_model.dl4j");
    	
    	try {
    		LibraryFile library=new LibraryFile();
    		library.openFile(new File("/Users/searleb/Documents/damien/hela_multiple_replicates_raws/2017aug23/23aug2017_hela_serum_timecourse_pool_wide_001_170829031834.dia.elib"));
    		ArrayList<LibraryEntry> entries=library.getAllEntries(false, aaConstants);

            MultiLayerNetwork model = ModelSerializer.restoreMultiLayerNetwork(saveLocation);
            System.out.println("Model loaded successfully.");
            
            long time=System.currentTimeMillis();
            
            ArrayList<INDArray> encodings=new ArrayList<INDArray>();
            
            for (LibraryEntry entry : entries) {
				encodings.add(encode(entry.getPeptideModSeq()));
			}
            
            INDArray inputData = Nd4j.concat(0, encodings.toArray(new INDArray[0]));
            inputData = inputData.reshape(entries.size(), MAXIMUM_PEPTIDE_LENGTH * AMINO_ACID_COUNT);
            INDArray output = model.output(inputData);

            for (int i = 0; i < entries.size(); i++) {
                System.out.println("Prediction for sequence " + entries.get(i).getPeptideModSeq() + "\t" + output.getDouble(i)+ "\t" + entries.get(i).getScanStartTime());
            }
            
    	} catch (Exception e) {
    		e.printStackTrace();
    	}
    }
    public static void main2(String[] args) {
    	File db=new File("/Users/searleb/Downloads/Chronologer_DB_220308.txt");
    	File saveLocation=new File(db.getParentFile(), "peptide_rt_model.dl4j");
    	
    	try {
    		LibraryFile library=new LibraryFile();
    		library.openFile(new File("/Users/searleb/Documents/damien/hela_multiple_replicates_raws/2017aug23/23aug2017_hela_serum_timecourse_pool_wide_001_170829031834.dia.elib"));
    		ArrayList<LibraryEntry> entries=library.getAllEntries(false, aaConstants);
    		
    	} catch (Exception e) {
    		e.printStackTrace();
    	}
    	
    	try {
            // Load the model
            MultiLayerNetwork model = ModelSerializer.restoreMultiLayerNetwork(saveLocation);
            System.out.println("Model loaded successfully.");
            
            long time=System.currentTimeMillis();
            for (int n = 0; n < 1000; n++) {	
	            ArrayList<INDArray> encodings=new ArrayList<INDArray>();
	            ArrayList<String> sequences=new ArrayList<String>();
	            
	            for (int i = 0; i < 1024; i++) {
					char[] aas=new char[RandomGenerator.randomIndex(22, i)+8];
					for (int j = 0; j < aas.length; j++) {
						aas[j]=AMINO_ACIDS.charAt(RandomGenerator.randomIndex(AMINO_ACID_COUNT, j));
					}
					String sequence=new String(aas);
					sequences.add(sequence);
					encodings.add(encode(sequence));			
				}
	            INDArray inputData = Nd4j.concat(0, encodings.toArray(new INDArray[0]));
	            inputData = inputData.reshape(sequences.size(), MAXIMUM_PEPTIDE_LENGTH * AMINO_ACID_COUNT);
	            INDArray output = model.output(inputData);
//	            for (int i = 0; i < sequences.size(); i++) {
//	                System.out.println("Prediction for sequence " + (i + 1) + ": " + output.getDouble(i));
//	            }
	            
	            System.out.println((n+1)+"--> "+(System.currentTimeMillis()-time)/1000f+" total sec, "+(System.currentTimeMillis()-time)/(n+1)+" msec per 1024");
			}
            
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void main3(String[] args) {
    	File db=new File("/Users/searleb/Downloads/Chronologer_DB_220308.txt");
    	File saveLocation=new File(db.getParentFile(), "peptide_rt_model.dl4j");
    	
    	ArrayList<ArrayList<DataSet>> dataSets = new ArrayList<>();
    	ArrayList<String> dataSetNames = new ArrayList<>();
    	
    	TableParser.parseTSV(db, new TableParserMuscle() {
        	ArrayList<DataSet> currentDataSet=null;
        	String currentDataSetName=null;
        	int count=0;
        	
			@Override
			public void processRow(Map<String, String> row) {
				String source=row.get("Source");
				if (!source.equals(currentDataSetName)) {
					
					currentDataSet=new ArrayList<DataSet>();
					dataSets.add(currentDataSet);
					dataSetNames.add(source);
					currentDataSetName=source;
				}
				
				String peptideModSeq=row.get("PeptideModSeq");
				float hirt=Float.parseFloat(row.get("HI"));

		        // One-hot encode the sequences
	            INDArray input = encode(peptideModSeq);
	            if (input==null) return;
	            
	            INDArray output = Nd4j.create(new double[][]{new double[] {hirt}});
	            currentDataSet.add(new DataSet(input, output));

				count++;
				if (count%10000==0) {
			    	System.out.println("Loading "+count+"...");
				}
			}
			
			@Override
			public void cleanup() {
			}
		});
    	System.out.println("Loaded "+dataSets.size()+" total datasets");

        MultiLayerNetwork model=createModel(MAXIMUM_PEPTIDE_LENGTH * AMINO_ACID_COUNT, 1, 
        		EMBED_DIMENSION, KERNEL_SIZE, N_RESNET_BLOCKS);
        System.out.println(model.summary());


    	@SuppressWarnings("unchecked")
		ArrayList<DataSet>[] trimmedDataSets = new ArrayList[dataSets.size()];
    	ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
        for(int e=0; e<NUM_EPOCHS; e++ ){
        	for (int d = 0; d < dataSets.size(); d++) {
            	System.out.println("Starting epoch "+e+" ("+dataSetNames.get(d)+")");
        		ArrayList<DataSet> dataSet;
        		if (e==0) {
        			dataSet=dataSets.get(d);
        		} else {
        			dataSet=trimmedDataSets[d];
        		}
        		
                DataSetIterator iterator = new ListDataSetIterator<>(dataSet, MAX_BATCH_SIZE);
                model.fit(iterator);

    			dataSet=dataSets.get(d);
                DataSetIterator fullIterator = new ListDataSetIterator<>(dataSet, MAX_BATCH_SIZE);
                INDArray results=model.output(fullIterator);
                
                TFloatArrayList deltas=new TFloatArrayList();
                for (int i = 0; i < dataSet.size(); i++) {
					double prediction=results.getDouble(i);
					double actual=dataSet.get(i).getLabels().getDouble(0);
					deltas.add((float)(actual-prediction));
				}
                float[] deltaArray=deltas.toArray();
                float low=QuickMedian.select(deltaArray, 0.005f);
                float high=QuickMedian.select(deltaArray, 0.995f);
                ArrayList<DataSet> trimmedData=new ArrayList<DataSet>();
                for (int i = 0; i < dataSet.size(); i++) {
					double prediction=results.getDouble(i);
					double actual=dataSet.get(i).getLabels().getDouble(0);
					double delta=actual-prediction;
					if (delta<high&&delta>low) {
						trimmedData.add(dataSet.get(i));
					}
				}
                trimmedDataSets[d]=trimmedData;
                
                ArrayList<XYPoint> histogram=PivotTableGenerator.createPivotTable(deltaArray, 25);
                traces.add(new XYTrace(histogram, GraphType.line, "Epoch "+e));
			}
        	
//        	if (e%5==0) {
//        		Charter.launchChart("Delta HI", "Count", true, traces.toArray(new XYTrace[0]));
//        	}
        }

        try {
        	ModelSerializer.writeModel(model, saveLocation, true);
        	System.out.println("Model saved successfully.");
        } catch (IOException ioe) {
        	ioe.printStackTrace();
        }
    }

    public static INDArray encode(String sequence) {
    	String[] aas=PeptideUtils.getAAs(sequence, aaConstants);
    	if (aas.length>MAXIMUM_PEPTIDE_LENGTH) return null;
    	
        INDArray encoded = Nd4j.zeros(MAXIMUM_PEPTIDE_LENGTH, AMINO_ACID_COUNT);
        
        for (int i = 0; i < aas.length; i++) {
            char aa = aas[i].charAt(0);
            if (aas[i].length()>1&&aa!='C') {
            	return null; // FIXME need to deal with other PTMs!
            }
            int index = AMINO_ACIDS.indexOf(aa);
            if (index >= 0) {
                encoded.putScalar(new int[]{i, index}, 1.0);
            } else {
                return null;
            }
        }
        return encoded.reshape(1, MAXIMUM_PEPTIDE_LENGTH * AMINO_ACID_COUNT);
    }
    
    public static MultiLayerNetwork createModel(int inputSize, int outputSize, int embedDim, int kernelSize, int resnetBlocks) {
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
        
        // dropout layer (keep 95% of changes)
        builder.layer(new DropoutLayer.Builder(0.95).build());
        
        // output layer
        builder.layer(new OutputLayer.Builder(LossFunctions.LossFunction.MSE)
                .nIn(embedDim)
                .nOut(outputSize)
                .activation(Activation.IDENTITY)
                .weightInit(WeightInit.XAVIER)
                .build());

        MultiLayerNetwork model = new MultiLayerNetwork(builder.build());
        model.init();
        model.setListeners(new ScoreIterationListener(10));

        return model;
    }
}