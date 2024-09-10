package edu.washington.gs.maccoss.encyclopedia.cli;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.StringTokenizer;
import java.util.TreeMap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaToPrositCSVParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaToPrositCSVParametersParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaFeaturePredictionModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.KoinaLibraryPredictionClient;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;

public class ConvertFastaToKoinaPrositLibrary {
	private static final String HTTPS_KOINA_WILHELMLAB_ORG_443 = "https://koina.wilhelmlab.org:443/";
	private static final String MODELS_FLAG = "-models";
	private static final String URL_FLAG = "-url";
	private static final String DELIM = ";";
	public static void main(String[] args) {
		HashMap<String, String> arguments= CommandLineParser.parseArguments(args);
		if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("CLI for Convert -> Create Prosit Library From FASTA using Koina");
			Logger.timelessLogLine("Required Parameters: ");
			Logger.timelessLogLine("\t-i\tinput .FASTA file");
			Logger.timelessLogLine("\nOther Parameters: ");
			Logger.timelessLogLine("\t-o\toutput .dlib file");

			TreeMap<String, String> defaults= new TreeMap<>(FastaToPrositCSVParametersParser.getDefaultParameters());
			int maxWidth=0;
			for (String key : defaults.keySet()) {
				if (key.length()>maxWidth) maxWidth=key.length();
			}
			maxWidth += 1;
			for (Entry<String, String> entry : defaults.entrySet()) {
				Logger.timelessLogLine("\t"+ General.formatCellToWidth(entry.getKey(), maxWidth)+" (default: "+entry.getValue()+")");
			}
			Logger.timelessLogLine("\t"+General.formatCellToWidth(MODELS_FLAG, maxWidth)+" (default: "+getModelParameter(KoinaFeaturePredictionModel.getDefaultModels())+")");
			Logger.timelessLogLine("\t"+ General.formatCellToWidth(URL_FLAG, maxWidth)+" (default: "+HTTPS_KOINA_WILHELMLAB_ORG_443+")");
			
			Logger.timelessLogLine("\nAvailable Models: ");

			Logger.timelessLogLine("\tFragmentation Models: ");
			for (KoinaFeaturePredictionModel model : KoinaFeaturePredictionModel.getFragmentationModels()) {
				Logger.timelessLogLine("\t\t"+model.getCodeName());
			}

			Logger.timelessLogLine("\tCCS Models: ");
			for (KoinaFeaturePredictionModel model : KoinaFeaturePredictionModel.getIMSModels()) {
				Logger.timelessLogLine("\t\t"+model.getCodeName());
			}

			Logger.timelessLogLine("\tRT Models: ");
			for (KoinaFeaturePredictionModel model : KoinaFeaturePredictionModel.getRTModels()) {
				Logger.timelessLogLine("\t\t"+model.getCodeName());
			}
			
		} else {
			convert(arguments);
		}
	}
	
	public static String getModelParameter(ArrayList<KoinaFeaturePredictionModel> models) {
		StringBuilder sb=new StringBuilder();
		for (KoinaFeaturePredictionModel model : models) {
			if (sb.length()>0) sb.append(DELIM);
			sb.append(model.getCodeName());
		}
		return sb.toString();
	}
	
	public static ArrayList<KoinaFeaturePredictionModel> getModels(String modelParameter) {
		if (modelParameter==null) {
			return KoinaFeaturePredictionModel.getDefaultModels();
		}
		
		StringTokenizer st=new StringTokenizer(modelParameter, DELIM);
		
		KoinaFeaturePredictionModel fragModel=null;
		KoinaFeaturePredictionModel imsModel=null;
		KoinaFeaturePredictionModel rtModel=null;
		
		while (st.hasMoreTokens()) {
			String modelName=st.nextToken();
			KoinaFeaturePredictionModel model=KoinaFeaturePredictionModel.getModel(modelName);
			if (KoinaFeaturePredictionModel.FRAGMENTATION_TYPE.equals(model.getModelType())) {
				fragModel=model;
			}
			if (KoinaFeaturePredictionModel.IMS_TYPE.equals(model.getModelType())) {
				imsModel=model;
			}
			if (KoinaFeaturePredictionModel.RT_TYPE.equals(model.getModelType())) {
				rtModel=model;
			}
		}
		
		if (fragModel==null) fragModel=KoinaFeaturePredictionModel.getDefaultModel(KoinaFeaturePredictionModel.FRAGMENTATION_TYPE);
		if (imsModel==null) fragModel=KoinaFeaturePredictionModel.getDefaultModel(KoinaFeaturePredictionModel.IMS_TYPE);
		if (rtModel==null) fragModel=KoinaFeaturePredictionModel.getDefaultModel(KoinaFeaturePredictionModel.RT_TYPE);
		
		ArrayList<KoinaFeaturePredictionModel> models=new ArrayList<KoinaFeaturePredictionModel>();
		models.add(fragModel);
		models.add(imsModel);
		models.add(rtModel);
		return models;
	}

	public static void convert(HashMap<String, String> arguments) {
		if (!arguments.containsKey("-i")) {
			Logger.errorLine("You are required to specify a FASTA file (-i), use -h for help");
			System.exit(1);
		}

		String modelNames = null;
		if (arguments.containsKey(MODELS_FLAG)) {
			modelNames = arguments.get(MODELS_FLAG);
		}
		ArrayList<KoinaFeaturePredictionModel> models = getModels(modelNames);

		String outputFile = null;
		if (arguments.containsKey("-o")) {
			outputFile = arguments.get("-o");
		}
		
		String baseURL=HTTPS_KOINA_WILHELMLAB_ORG_443;
		if (arguments.containsKey(URL_FLAG)) {
			baseURL = arguments.get(URL_FLAG);
		}

		File fastaFile = new File(arguments.get("-i"));
		FastaToPrositCSVParameters params= FastaToPrositCSVParametersParser.parseParameters(arguments);
		try {
			if (fastaFile.exists()) {
				KoinaLibraryPredictionClient.writeLibrary(baseURL, models, outputFile, fastaFile, params.getEnzyme(), params.getDefaultNCE(), params.getDefaultCharge(),
						params.getMinCharge(), params.getMaxCharge(), params.getMaxMissedCleavage(), new Range(params.getMinMz(), params.getMaxMz()), 
						params.isAdjustNCEForDIA(), params.isAddDecoys(), SearchParameterParser.getDefaultParametersObject(), new EmptyProgressIndicator(false));
				
				Logger.logLine("Finished conversion from FASTA to Prosit dlib");
			} else {
				Logger.logLine("You must specify a FASTA file!");
			}
		} catch (Exception e) {
			Logger.errorLine("Encountered Fatal Error!");
			Logger.errorException(e);
		}
	}
}
