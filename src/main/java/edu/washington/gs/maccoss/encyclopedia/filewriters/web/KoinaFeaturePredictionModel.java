package edu.washington.gs.maccoss.encyclopedia.filewriters.web;

import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import edu.washington.gs.maccoss.encyclopedia.algorithms.prediction.AminoAcidEncoding;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation.Prosit2019HCDModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation.Prosit2020CIDModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation.Prosit2020HCDModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation.Prosit2020TMTModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.fragmentation.Prosit2023timsTOFModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.ims.AlphaPeptDeepIMSModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.ims.IM2DeepIMSModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt.AlphaPeptDeepRTModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt.DeepLCHelaRTModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt.Prosit2019RTModel;
import edu.washington.gs.maccoss.encyclopedia.filewriters.web.models.rt.Prosit2020TMTRTModel;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public interface KoinaFeaturePredictionModel {
	public static final String FRAGMENTATION_TYPE="FRAG";
	public static final String IMS_TYPE="CCS";
	public static final String RT_TYPE="RT";
	
	public String getName();
	public String getCodeName();
	public String getModelType();
	public URL getURL();
	public void updatePeptides(List<KoinaPrecursor> peptides);
	public boolean canModelPeptide(AminoAcidEncoding[] aas, byte precursorCharge);

	public static ArrayList<KoinaFeaturePredictionModel> getDefaultModels() {
		ArrayList<KoinaFeaturePredictionModel> models=new ArrayList<KoinaFeaturePredictionModel>();
		models.add(getDefaultModel(FRAGMENTATION_TYPE));
		models.add(getDefaultModel(IMS_TYPE));
		models.add(getDefaultModel(RT_TYPE));
		return models;
	}

	public static KoinaFeaturePredictionModel getDefaultModel(String modelType) {
		if (FRAGMENTATION_TYPE.equals(modelType)) {
			return new Prosit2020HCDModel();
		}
		if (IMS_TYPE.equals(modelType)) {
			return new IM2DeepIMSModel();
		}
		if (RT_TYPE.equals(modelType)) {
			return new Prosit2019RTModel();
		}

		throw new EncyclopediaException("Unexpected type ["+modelType+"]");
	}
	
	public static KoinaFeaturePredictionModel getModel(String modelName, String modelType) {
		if (FRAGMENTATION_TYPE.equals(modelType)) {
			for (KoinaFeaturePredictionModel model : getFragmentationModels()) {
				if (model.getCodeName().equalsIgnoreCase(modelName)) return model;
			}
		}
		if (IMS_TYPE.equals(modelType)) {
			for (KoinaFeaturePredictionModel model : getIMSModels()) {
				if (model.getCodeName().equalsIgnoreCase(modelName)) return model;
			}
		}
		if (RT_TYPE.equals(modelType)) {
			for (KoinaFeaturePredictionModel model : getRTModels()) {
				if (model.getCodeName().equalsIgnoreCase(modelName)) return model;
			}
		}
		
		throw new EncyclopediaException("Unexpected model name: ["+modelName+"] for type ["+modelType+"]");
	}
	
	public static KoinaFeaturePredictionModel getModel(String modelName) {
		for (KoinaFeaturePredictionModel model : getFragmentationModels()) {
			if (model.getCodeName().equalsIgnoreCase(modelName)) return model;
		}
		for (KoinaFeaturePredictionModel model : getIMSModels()) {
			if (model.getCodeName().equalsIgnoreCase(modelName)) return model;
		}
		for (KoinaFeaturePredictionModel model : getRTModels()) {
			if (model.getCodeName().equalsIgnoreCase(modelName)) return model;
		}
		
		throw new EncyclopediaException("Unexpected model name: ["+modelName+"]");
	}
	
	/**
	 * default model is first
	 * @return
	 */
	public static ArrayList<KoinaFeaturePredictionModel> getFragmentationModels() {
		ArrayList<KoinaFeaturePredictionModel> models=new ArrayList<KoinaFeaturePredictionModel>();
		models.add(new Prosit2020HCDModel());
		models.add(new Prosit2020CIDModel());
		models.add(new Prosit2023timsTOFModel());
		models.add(new Prosit2019HCDModel());
		models.add(new Prosit2020TMTModel(true));
		models.add(new Prosit2020TMTModel(false));
		return models;
	}
	
	/**
	 * default model is first
	 * @return
	 */
	public static ArrayList<KoinaFeaturePredictionModel> getIMSModels() {
		ArrayList<KoinaFeaturePredictionModel> models=new ArrayList<KoinaFeaturePredictionModel>();
		models.add(new IM2DeepIMSModel());
		models.add(new AlphaPeptDeepIMSModel());
		return models;
	}
	
	/**
	 * default model is first
	 * @return
	 */
	public static ArrayList<KoinaFeaturePredictionModel> getRTModels() {
		ArrayList<KoinaFeaturePredictionModel> models=new ArrayList<KoinaFeaturePredictionModel>();
		models.add(new Prosit2019RTModel());
		models.add(new DeepLCHelaRTModel());
		models.add(new AlphaPeptDeepRTModel());
		models.add(new Prosit2020TMTRTModel());
		return models;
	}
}
