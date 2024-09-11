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
	public URL getURL(String baseURL);
	public void updatePeptides(List<KoinaPrecursor> peptides, String baseURL);
	public boolean canModelPeptide(AminoAcidEncoding[] aas, byte precursorCharge);

	/**
	 * Retrieves the default set of Koina feature prediction models. This includes a fragmentation model,
	 * an ion mobility spectrometry (IMS) model, and a retention time (RT) model.
	 *
	 * @return an ArrayList of the default KoinaFeaturePredictionModel objects.
	 */
	public static ArrayList<KoinaFeaturePredictionModel> getDefaultModels() {
		ArrayList<KoinaFeaturePredictionModel> models=new ArrayList<KoinaFeaturePredictionModel>();
		models.add(getDefaultModel(FRAGMENTATION_TYPE));
		models.add(getDefaultModel(IMS_TYPE));
		models.add(getDefaultModel(RT_TYPE));
		return models;
	}

	/**
	 * Retrieves the default model for a specific model type (e.g., FRAG, CCS, or RT).
	 *
	 * @param modelType the type of the model (e.g., FRAGMENTATION_TYPE, IMS_TYPE, or RT_TYPE).
	 * @return the default KoinaFeaturePredictionModel corresponding to the modelType.
	 * @throws EncyclopediaException if an unexpected model type is provided.
	 */
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
	
	/**
	 * Retrieves a specific KoinaFeaturePredictionModel by its name and type.
	 *
	 * @param modelName the name of the model to retrieve.
	 * @param modelType the type of the model (e.g., FRAGMENTATION_TYPE, IMS_TYPE, or RT_TYPE).
	 * @return the KoinaFeaturePredictionModel corresponding to the modelName and modelType.
	 * @throws EncyclopediaException if the model name or type is unexpected or not found.
	 */
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
	
	/**
	 * Retrieves a specific KoinaFeaturePredictionModel by its name. This method searches across all
	 * types of models (fragmentation, IMS, and RT).
	 *
	 * @param modelName the name of the model to retrieve.
	 * @return the KoinaFeaturePredictionModel corresponding to the modelName.
	 * @throws EncyclopediaException if the model name is unexpected or not found.
	 */
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
	 * Retrieves a list of all available fragmentation models, with the default model listed first.
	 *
	 * @return an ArrayList of KoinaFeaturePredictionModel objects for fragmentation.
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
	 * Retrieves a list of all available IMS (ion mobility spectrometry) models, with the default model listed first.
	 *
	 * @return an ArrayList of KoinaFeaturePredictionModel objects for IMS.
	 */
	public static ArrayList<KoinaFeaturePredictionModel> getIMSModels() {
		ArrayList<KoinaFeaturePredictionModel> models=new ArrayList<KoinaFeaturePredictionModel>();
		models.add(new IM2DeepIMSModel());
		models.add(new AlphaPeptDeepIMSModel());
		return models;
	}
	
	/**
	 * Retrieves a list of all available retention time (RT) models, with the default model listed first.
	 *
	 * @return an ArrayList of KoinaFeaturePredictionModel objects for RT.
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
