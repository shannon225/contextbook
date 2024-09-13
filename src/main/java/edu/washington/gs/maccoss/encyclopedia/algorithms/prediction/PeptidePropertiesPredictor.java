package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;


import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.zip.DataFormatException;

import org.deeplearning4j.nn.conf.ComputationGraphConfiguration;
import org.deeplearning4j.nn.conf.ConvolutionMode;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.graph.ElementWiseVertex;
import org.deeplearning4j.nn.conf.graph.MergeVertex;
import org.deeplearning4j.nn.conf.graph.ReshapeVertex;
import org.deeplearning4j.nn.conf.inputs.InputType;
import org.deeplearning4j.nn.conf.layers.ActivationLayer;
import org.deeplearning4j.nn.conf.layers.BatchNormalization;
import org.deeplearning4j.nn.conf.layers.ConvolutionLayer;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.GlobalPoolingLayer;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
import org.deeplearning4j.nn.conf.layers.PoolingType;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.nn.weights.WeightInit;
import org.deeplearning4j.optimize.listeners.ScoreIterationListener;
import org.deeplearning4j.util.ModelSerializer;
import org.nd4j.linalg.activations.Activation;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.nd4j.linalg.learning.config.Adam;
import org.nd4j.linalg.lossfunctions.LossFunctions;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Correlation;
import edu.washington.gs.maccoss.encyclopedia.utils.math.PivotTableGenerator;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class PeptidePropertiesPredictor {
	private static final int REPORTING_BATCH_SIZE = 1000;
	private static final int INITIAL_BATCH_SIZE = 64;
	public static final int MAX_BATCH_SIZE = 1024;
	private static final int EPOCHS_TO_2X_BATCH = 99999; // never increase
    public static final int EMBED_DIMENSION = 128;
    public static final int N_RESNET_BLOCKS = 3;
    public static final int KERNEL_SIZE = 9;
    public static final int NUM_EPOCHS = 3; // only limited epochs with iterative training
    
    public static void main(String[] args) throws Exception {
    	File dir=new File("C:\\Users\\searl\\Documents\\projects\\prosit_examples_final\\");
    	
    	// FIXME
    	dir=new File("/Users/searleb/Downloads/");
    	File saveLocation=new File(dir, "peptide_prediction_model.dl4j");
    	if (saveLocation.exists()) {
    		saveLocation.delete();
    	}
    	
    	String[] datasetOrder=new String[] {
        		"UP000005640_9606.fasta.trypsin.z2_nce33.dlib.z3_nce33.dlib",
        		"UP000005640_9606.fasta.trypsin.z3_nce33.dlib.z3_nce33.dlib",
        		"UP000005640_9606.fasta.trypsin.z4_nce33.dlib.z3_nce33.dlib",
        		"UP000005640_9606.fasta.trypsin.z1_nce33.dlib.z3_nce33.dlib",
        		"UP000005640_9606.fasta.trypsin.z5_nce33.dlib.z3_nce33.dlib",
        		"UP000005640_9606.fasta.trypsin.z6_nce33.dlib.z3_nce33.dlib",
        		"UP000005640_9606.fasta.lys-n.z3_nce33.dlib.z3_nce33.dlib",
        		"UP000005640_9606.fasta.asp-n.z3_nce33.dlib.z3_nce33.dlib",
        		"UP000005640_9606.fasta.chymotrypsin.z3_nce33.dlib.z3_nce33.dlib",
        		"UP000005640_9606.fasta.glu-c.z3_nce33.dlib.z3_nce33.dlib",
        };
    	
    	// FIXME
    	datasetOrder=new String[] {"trimmed_Trypsin.csv.dlib"};
    	
    	try {
    		// start on the "easiest" libraries first	
    		for (int i = 0; i < datasetOrder.length; i++) {
				File[] listFiles=new File[i+1];
				for (int j = 0; j < listFiles.length; j++) {
					listFiles[j]=new File(dir, datasetOrder[j]);
				}

	    		runCycle(listFiles, saveLocation);
			}
    		
    	} catch (IOException | SQLException | DataFormatException e) {
    		Logger.errorLine("Found unexpected IO exception!");
    		Logger.errorException(e);
    	}
    }
    
    private static void runCycle(File[] listFiles, File saveLocation) throws IOException, SQLException, DataFormatException {
    	
    	SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
    	AminoAcidConstants aaConstants=parameters.getAAConstants();

		HashMap<String, PeptideEncoding> peptideEncodingBySequence=new HashMap<String, PeptideEncoding>();
    	
		for (File f : listFiles) {
        	TObjectFloatHashMap<String> rtByPeptideModSeq=new TObjectFloatHashMap<String>();
    		String rtName=f.getName().substring(0, f.getName().length()-".z3_nce33.dlib".length())+".txt_rts.txt";
    		File rtFile = new File(f.getParent(), rtName);
    		if (rtFile.exists()) {
	    		Logger.logLine("Reading RT file: "+rtName);
	
				TableParser.parseTSV(rtFile, new TableParserMuscle() {
					
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
    		
    		Logger.logLine("Reading Library file: "+f.getName());
    		LibraryFile library=new LibraryFile();
    		library.openFile(f);
    		ArrayList<LibraryEntry> entries=library.getAllEntries(false, aaConstants);
    		for (LibraryEntry entry : entries) {
    			if (entry.getPrecursorCharge()>PeptideEncoding.MAX_CHARGE||entry.getPeptideSeq().length()>(PeptideEncoding.MAX_PEPTIDE_LENGTH-2)) {
    				continue;
    			}
				String key=entry.getPeptideModSeq()+"_"+entry.getPrecursorCharge();
				if (rtByPeptideModSeq.contains(entry.getPeptideModSeq())) {
					PeptideEncoding encoding=new PeptideEncoding(entry, rtByPeptideModSeq.get(entry.getPeptideModSeq()), parameters);
					peptideEncodingBySequence.put(key, encoding);
				} else {
					PeptideEncoding encoding=new PeptideEncoding(entry, entry.getRetentionTimeInSec(), parameters);
					peptideEncodingBySequence.put(key, encoding);
				}
			}
    		library.close();
    	}

		Runtime instance=Runtime.getRuntime();
		long usedMemory=instance.totalMemory()-instance.freeMemory();
	 	Logger.logLine((usedMemory/1048576)+" of "+(instance.maxMemory()/1048576)+" MB used");

    	Logger.logLine("Loaded data "+peptideEncodingBySequence.size()+" total precursors.");

    	ComputationGraph model;
        if (saveLocation.exists()&&saveLocation.canRead()) {
            model = ModelSerializer.restoreComputationGraph(saveLocation);
            model.getConfiguration().setIterationCount(0);
            model.getConfiguration().setEpochCount(0);
            model.setListeners(new ScoreIterationListener(REPORTING_BATCH_SIZE));
            Logger.logLine("Existing model loaded successfully.");
        } else {
        	// create new model
	        model=createModel(EMBED_DIMENSION, KERNEL_SIZE, N_RESNET_BLOCKS, 0.9f);
	        //model=ResNetBuilder.createModel(EMBED_DIMENSION, KERNEL_SIZE, N_RESNET_BLOCKS, 0.9f);
            model.setListeners(new ScoreIterationListener(REPORTING_BATCH_SIZE));
            Logger.logLine("New model created.");
        }
        Logger.logLine(model.summary());
        
        int batchSize=INITIAL_BATCH_SIZE;
    	double bestAverageContrastAngle=-Double.MAX_VALUE;
    	
    	ArrayList<PeptideEncoding> dataSet=new ArrayList<PeptideEncoding>(peptideEncodingBySequence.values());
    	ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
        for(int e=1; e<=NUM_EPOCHS; e++ ){
        	if (e%EPOCHS_TO_2X_BATCH==0) {
        		batchSize=batchSize*2;
        	}
        	Logger.logLine("Starting epoch "+e+", batch size: "+batchSize);
        	Collections.shuffle(dataSet);

        	MultiDataSetIterator iterator = new EncodedMultiDataSetIterator(dataSet, batchSize);
            model.fit(iterator);
            
            MultiDataSetIterator fullIterator = new EncodedMultiDataSetIterator(dataSet, MAX_BATCH_SIZE);
            INDArray[] results=model.output(fullIterator);
            
            INDArray fragmentationResults=results[0].reshape(dataSet.size(), PeptideEncoding.ENCODED_OUTPUT_FRAGMENT_SIZE);
            INDArray rtResults=results[1].reshape(dataSet.size(), PeptideEncoding.ENCODED_OUTPUT_RT_SIZE);
            INDArray imsResults=results[2].reshape(dataSet.size(), PeptideEncoding.ENCODED_OUTPUT_IMS_SIZE);

            TFloatArrayList contrastAngles=new TFloatArrayList();
            for (int i = 0; i < dataSet.size(); i++) {
    			INDArray fragRow=fragmentationResults.getRow(i);
    			INDArray rtRow=rtResults.getRow(i);
    			INDArray imsRow=imsResults.getRow(i);
				PeptideEncoding peptideEncoding = dataSet.get(i);
				//float spectralContrastAngle=peptideEncoding.score(row);
				LibraryEntry predicted=peptideEncoding.outputToEntry(new HashSet<>(), new INDArray[] {fragRow, rtRow, imsRow}, aaConstants);
				LibraryEntry original=peptideEncoding.encodingToEntry(new HashSet<>(), aaConstants);
				contrastAngles.add((float)Correlation.getSpectralAngle(predicted, original, parameters.getFragmentTolerance()));
			}
        	double averageContrastAngle=0.0;
        	for (int i = 0; i < contrastAngles.size(); i++) {
				averageContrastAngle+=contrastAngles.get(i);
			}
        	averageContrastAngle=averageContrastAngle/contrastAngles.size();
            
        	Logger.logLine("Epoch "+e+" mean contrast angle: "+averageContrastAngle+" from "+contrastAngles.size()+" total peptides");

        	if (averageContrastAngle>bestAverageContrastAngle) {
        		bestAverageContrastAngle=averageContrastAngle;
	            try {
	            	ModelSerializer.writeModel(model, saveLocation, true);
	            	Logger.logLine("Model saved successfully.");
	            } catch (IOException ioe) {
	            	ioe.printStackTrace();
	            }
        	} else {
            	Logger.logLine("Final model wasn't better, so skipping temporary save.");
        	}
        	
        	XYTrace thisTrace=new XYTrace(PivotTableGenerator.createPivotTable(contrastAngles.toArray(), -1f, 1f, 0.01f), GraphType.area, "Epoch "+e);
        	traces.add(thisTrace);
        	
        	Charter.launchChart(Charter.getChart("Epoch", "ContrastAngle", true, traces.toArray(new XYTrace[0])), "Traces");
        }
	}

	private static ComputationGraph createModel(int embedDim, int kernelSize, int resnetBlocks, double dropRate) {
		embedDim=32;
		ComputationGraphConfiguration.GraphBuilder graph = new NeuralNetConfiguration.Builder()
				.weightInit(WeightInit.XAVIER).updater(new Adam(0.001)).graphBuilder().addInputs("sequence", "charge")
				.setInputTypes(
						InputType.feedForward(
								AminoAcidEncoding.MAX_ENCODING_LENGTH * PeptideEncoding.MAX_PEPTIDE_LENGTH),
						InputType.feedForward(PeptideEncoding.ENCODED_INPUT_CHARGE_SIZE));

		// Sequence embedding
		graph.addLayer("seq_embedding", new DenseLayer.Builder()
				.nIn(AminoAcidEncoding.MAX_ENCODING_LENGTH * PeptideEncoding.MAX_PEPTIDE_LENGTH).nOut(embedDim).build(),
				"sequence");

		// Charge embedding
		graph.addLayer("charge_embedding", new DenseLayer.Builder().nIn(PeptideEncoding.MAX_CHARGE).nOut(embedDim)
				.activation(Activation.IDENTITY).build(), "charge");

		// Combine embeddings
		graph.addVertex("merged", new MergeVertex(), "seq_embedding", "charge_embedding");

		// Reshape merged output to [batchSize, channels=1, height=1, width=2*embedDim]
		graph.addVertex("reshaped", new ReshapeVertex(-1, 1, 1, 2 * embedDim), "merged");

		String previousLayer = "reshaped";
		int[] filters = { 32, 32, 64 };
		int[] kernelSizes = { kernelSize, 1 };
		int[] strides = { 1, 1 };

		for (int i = 1; i <= resnetBlocks; i++) {
			previousLayer = convBlock(graph, kernelSizes, filters, Integer.toString(i), strides, previousLayer,
					dropRate);
		}

		graph.addLayer("global_pooling", new GlobalPoolingLayer.Builder().poolingType(PoolingType.AVG).build(),
				previousLayer);
		previousLayer = "global_pooling";

		int numChannelsAfterPooling = filters[2];

		// Output layers
		graph.addLayer("fragmentation_output",
				new OutputLayer.Builder(LossFunctions.LossFunction.MSE).activation(Activation.IDENTITY)
						.nIn(numChannelsAfterPooling).nOut(PeptideEncoding.ENCODED_OUTPUT_FRAGMENT_SIZE).build(),
				previousLayer);

		// Retention time output branch
		graph.addLayer("ret_time_output",
				new OutputLayer.Builder(LossFunctions.LossFunction.MSE).activation(Activation.IDENTITY)
						.nIn(numChannelsAfterPooling).nOut(PeptideEncoding.ENCODED_OUTPUT_RT_SIZE).build(),
				previousLayer);

		// Ion mobility output branch
		graph.addLayer("ion_mobility_output",
				new OutputLayer.Builder(LossFunctions.LossFunction.MSE).activation(Activation.IDENTITY)
						.nIn(numChannelsAfterPooling).nOut(PeptideEncoding.ENCODED_OUTPUT_IMS_SIZE).build(),
				previousLayer);

		graph.setOutputs("fragmentation_output", "ret_time_output", "ion_mobility_output");

		ComputationGraph model = new ComputationGraph(graph.build());
		model.init();
		return model;
	}

	public static String convBlock(ComputationGraphConfiguration.GraphBuilder graph, int[] kernelSize, int[] filters,
			String stage, int[] stride, String input, double dropRate) {
		String convName = "res" + stage + "_branch";
		String batchName = "bn" + stage + "_branch";
		String activationName = "act" + stage + "_branch";
		String shortcutName = "short" + stage + "_branch";

		// Convolutional layers with dropout
		graph.addLayer(convName + "2a",
				new ConvolutionLayer.Builder(kernelSize, stride).nOut(filters[0]).convolutionMode(ConvolutionMode.Same)
						.dropOut(dropRate).build(),
				input)
				.addLayer(batchName + "2a", new BatchNormalization.Builder().nOut(filters[0]).build(), convName + "2a")
				.addLayer(activationName + "2a", new ActivationLayer.Builder().activation(Activation.RELU).build(),
						batchName + "2a")

				.addLayer(convName + "2b",
						new ConvolutionLayer.Builder(kernelSize, new int[] { 1, 1 }).nOut(filters[1])
								.convolutionMode(ConvolutionMode.Same).dropOut(dropRate).build(),
						activationName + "2a")
				.addLayer(batchName + "2b", new BatchNormalization.Builder().nOut(filters[1]).build(), convName + "2b")
				.addLayer(activationName + "2b", new ActivationLayer.Builder().activation(Activation.RELU).build(),
						batchName + "2b")

				.addLayer(convName + "2c",
						new ConvolutionLayer.Builder(new int[] { 1, 1 }, new int[] { 1, 1 }).nOut(filters[2])
								.convolutionMode(ConvolutionMode.Same).dropOut(dropRate).build(),
						activationName + "2b")
				.addLayer(batchName + "2c", new BatchNormalization.Builder().nOut(filters[2]).build(), convName + "2c")

				// Shortcut connection
				.addLayer(convName + "1",
						new ConvolutionLayer.Builder(new int[] { 1, 1 }, stride).nOut(filters[2])
								.convolutionMode(ConvolutionMode.Same).dropOut(dropRate).build(),
						input)
				.addLayer(batchName + "1", new BatchNormalization.Builder().nOut(filters[2]).build(), convName + "1")

				.addVertex(shortcutName, new ElementWiseVertex(ElementWiseVertex.Op.Add), batchName + "2c",
						batchName + "1")
				.addLayer(convName, new ActivationLayer.Builder().activation(Activation.RELU).build(), shortcutName);

		return convName;
	}
}