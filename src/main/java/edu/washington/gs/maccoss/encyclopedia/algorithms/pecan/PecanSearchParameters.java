package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import edu.washington.gs.maccoss.encyclopedia.Pecanpie;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorVersion;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ScoringBreadthType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.parameters.InstrumentSpecificSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

//@Immutable
public class PecanSearchParameters extends SearchParameters {
	public static final boolean DEFAULT_HANDLE_N_TERM_METHIONINE = true; // TODO: allow user to set this parameter
	private final int minPeptideLength;
	private final int maxPeptideLength;
	private final int maxMissedCleavages;
	private final byte minCharge;
	private final byte maxCharge;
	private final int numberOfReportedPeaks;
	private final boolean addDecoysToBackgound;
	private final boolean dontRunDecoys; // only for testing
	private final float alpha;
	private final float beta;
	private final boolean requireVariableMods;
	
	public String toString() {
		final StringBuilder sb=new StringBuilder(super.toString());
		sb.append(" -minLength "+minPeptideLength+"\n");
		sb.append(" -maxLength "+maxPeptideLength+"\n");
		sb.append(" -maxMissedCleavage "+maxMissedCleavages+"\n");
		sb.append(" -minCharge "+minCharge+"\n");
		sb.append(" -maxCharge "+maxCharge+"\n");
		sb.append(" -numberOfReportedPeaks "+numberOfReportedPeaks+"\n");
		sb.append(" -addDecoysToBackground "+addDecoysToBackgound+"\n");
		sb.append(" -dontRunDecoys "+dontRunDecoys+"\n");
		sb.append(" -alpha "+alpha+"\n");
		sb.append(" -beta "+beta+"\n");
		sb.append(" -requireVariableMods "+requireVariableMods+"\n");
		
		return sb.toString();
	}
	
	@Override
	public HashMap<String, String> toParameterMap() {
		HashMap<String, String> map=super.toParameterMap();
		map.put("-minLength", minPeptideLength+"");
		map.put("-maxLength", maxPeptideLength+"");
		map.put("-maxMissedCleavage", maxMissedCleavages+"");
		map.put("-minCharge", minCharge+"");
		map.put("-maxCharge", maxCharge+"");
		map.put("-expectedPeakWidth", expectedPeakWidth+"");
		map.put("-numberOfReportedPeaks", numberOfReportedPeaks+"");
		map.put("-addDecoysToBackground", addDecoysToBackgound+"");
		map.put("-dontRunDecoys", dontRunDecoys+"");
		map.put("-alpha", alpha+"");
		map.put("-beta", beta+"");
		map.put("-requireVariableMods", "false");
		return map;
	}
	
	@Override
	public HashMap<String, String> getNonDefaultParameters() {
		HashMap<String, String> map=toParameterMap();
		HashMap<String, String> defaults=PecanParameterParser.getDefaultParametersObject().toParameterMap();
		HashMap<String, String> nonDefaults=new HashMap<>();
		
		for (Entry<String, String> entry : map.entrySet()) {
			if (defaults.containsKey(entry.getKey())) {
				String defaultValue=defaults.get(entry.getKey());
				if (!entry.getValue().equals(defaultValue)) {
					nonDefaults.put(entry.getKey(), entry.getValue());
				}
			} else {
				nonDefaults.put(entry.getKey(), entry.getValue());
			}
		}
		return nonDefaults;
	}
	
	@Override
	public void writeToXML(Document doc, Element parentElement) {
		HashMap<String, String> nondefaults=getNonDefaultParameters();
		
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);
		
		for (Entry<String, String> entry : nondefaults.entrySet()) {
			Element param=doc.createElement("param");
			param.setAttribute("key", entry.getKey());
			param.setAttribute("value", entry.getValue());
			rootElement.appendChild(param);
		}
	}
	
	public static PecanSearchParameters readFromXML(Document doc, Element rootElement) {
		if (!rootElement.getTagName().equals(PecanSearchParameters.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+PecanSearchParameters.class.getSimpleName()+"]");
		}
		
		HashMap<String, String> paramMap=PecanParameterParser.getDefaultParametersObject().toParameterMap();
		
		NodeList nodes=rootElement.getChildNodes();
		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
            if (node.getNodeType() == Node.ELEMENT_NODE) {
                Element element = (Element) node;
                if ("param".equals(element.getTagName())) {
                	String key=element.getAttribute("key");
                	String value=element.getAttribute("value");
                	paramMap.put(key, value);
                }
            }
		}
		return PecanParameterParser.parseParameters(paramMap);
	}
	
	@Override
	public void savePreferences(File backgroundFastaFile, File targetFastaFile) throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("pecan");
		HashMap<String, String> map=toParameterMap();
		if (backgroundFastaFile!=null) map.put(Pecanpie.BACKGROUND_FASTA_TAG, backgroundFastaFile.getAbsolutePath());
		if (targetFastaFile!=null) map.put(Pecanpie.TARGET_FASTA_TAG, targetFastaFile.getAbsolutePath());
		for (Entry<String, String> entry : map.entrySet()) {
			//System.out.println("Writing Pecan preference "+entry.getKey()+" = "+entry.getValue());
			prefs.put(entry.getKey(), entry.getValue());
		}
		prefs.flush();
	}
	
	public static HashMap<String, String> readPreferences() throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("pecan");
		HashMap<String, String> map=new HashMap<String, String>();
		for (String key : prefs.keys()) {
			String value=prefs.get(key, "");
			//System.out.println("Reading Pecan preference "+key+" = "+value);
			map.put(key, value);
		}
		return map;
	}
	
	/** used by CLI
	 */
	public PecanSearchParameters(
			AminoAcidConstants aaConstants, 
			FragmentationType fragType, 
			MassTolerance precursorTolerance, 
			double precursorOffsetPPM, 
			double precursorIsolationMargin, 
			MassTolerance fragmentTolerance,
			double fragmentOffsetPPM, 
			DigestionEnzyme enzyme, 
			float expectedPeakWidth,
			int minPeptideLength, 
			int maxPeptideLength, 
			int maxMissedCleavages, 
			byte minCharge, 
			byte maxCharge, 
			int numberOfReportedPeaks,
			boolean addDecoysToBackgound,
			boolean dontRunDecoys,
			float percolatorThreshold,
			float percolatorProteinThreshold,
			boolean usePercolator,
			PercolatorVersion percolatorVersionNumber,
			int percolatorTrainingSetSize,
			float percolatorTrainingSetThreshold,
			int percolatorTrainingIterations,
			float alpha,
			float beta,
			DataAcquisitionType dataAcquisitionType,
			int numberOfThreadsUsed,
			float targetWindowCenter,
			float precursorWindowSize, 
			int numberOfQuantitativePeaks, 
			int minNumOfQuantitativePeaks, 
			int topNTargetsUsed,
			float minIntensity, 
			float minIntensityNumIons, 
			boolean quantifyAcrossSamples, 
			boolean verifyModificationIons, 
			boolean requireVariableMods, 
			float rtWindowInMin,
			int minNumIntegratedRTPoints,
			boolean filterPeaklists, 
			boolean doNotUseGlobalFDR, 
			Optional<File> precursorIsolationRangeFile, 
			Optional<File> percolatorModelFile, 
			boolean normalizeByTIC,
			boolean subtractBackground,
			boolean maskBadIntegrations,
			boolean adjustInferredRTBoundaries,
			boolean integratePrecursors,
			InstrumentSpecificSearchParameters instrument,
			boolean enableAdvancedOptions) {
		super(
				aaConstants,
				fragType,
				precursorTolerance,
				precursorOffsetPPM,
				precursorIsolationMargin,
				fragmentTolerance,
				fragmentOffsetPPM,
				fragmentTolerance,
				enzyme,
				percolatorThreshold,
				percolatorProteinThreshold,
				usePercolator,
				percolatorVersionNumber,
				percolatorTrainingSetSize,
				percolatorTrainingSetThreshold,
				percolatorTrainingIterations,
				dataAcquisitionType,
				numberOfThreadsUsed,
				expectedPeakWidth,
				targetWindowCenter,
				precursorWindowSize,
				numberOfQuantitativePeaks,
				minNumOfQuantitativePeaks,
				topNTargetsUsed,
				minIntensity,
				minIntensityNumIons,
				Optional.ofNullable((PeptideModification)null),
				ScoringBreadthType.ENTIRE_RT_WINDOW,
				0,
				quantifyAcrossSamples,
				verifyModificationIons,
				rtWindowInMin,
				minNumIntegratedRTPoints,
				filterPeaklists,
				doNotUseGlobalFDR, 
				precursorIsolationRangeFile,
				percolatorModelFile,
				normalizeByTIC,
				subtractBackground,
				maskBadIntegrations,
				adjustInferredRTBoundaries,
				integratePrecursors,
				instrument,
				enableAdvancedOptions
		);
		this.minPeptideLength=minPeptideLength;
		this.maxPeptideLength=maxPeptideLength;
		this.maxMissedCleavages=maxMissedCleavages;
		this.minCharge=minCharge;
		this.maxCharge=maxCharge;
		this.numberOfReportedPeaks=numberOfReportedPeaks;
		this.addDecoysToBackgound=addDecoysToBackgound;
		this.dontRunDecoys=dontRunDecoys;
		this.alpha=alpha;
		this.beta=beta;
		this.requireVariableMods=requireVariableMods;
	}

	/** used by GUI
	 */
	public PecanSearchParameters(
			AminoAcidConstants aaConstants,
			FragmentationType fragType,
			MassTolerance precursorTolerance,
			MassTolerance fragmentTolerance,
			DigestionEnzyme enzyme,
			PercolatorVersion percolatorVersionNumber,
			float percolatorThreshold,
			float percolatorProteinThreshold,
			int percolatorTrainingSetSize,
			float percolatorTrainingSetThreshold,
			int percolatorTrainingIterations,
			int maxMissedCleavages,
			byte minCharge,
			byte maxCharge,
			DataAcquisitionType dataAcquisitionType,
			float precursorWindowSize,
			int numberOfJobs,
			int numberOfQuantitativePeaks,
			int minNumOfQuantitativePeaks,
			int topNTargetsUsed,
			float minIntensity,
			float numberOfExtraDecoyLibrariesSearched,
			boolean quantifyAcrossSamples,
			boolean verifyModificationIons,
			boolean requireVariableMods,
			InstrumentSpecificSearchParameters instrument
	) {
		super(
				aaConstants,
				fragType,
				precursorTolerance,
				0.0,
				0.0,
				fragmentTolerance,
				0.0,
				fragmentTolerance,
				enzyme,
				percolatorThreshold,
				percolatorProteinThreshold,
				true,
				percolatorVersionNumber,
				percolatorTrainingSetSize,
				percolatorTrainingSetThreshold,
				percolatorTrainingIterations,
				dataAcquisitionType,
				numberOfJobs,
				24f,
				-1f,
				precursorWindowSize,
				numberOfQuantitativePeaks,
				minNumOfQuantitativePeaks,
				topNTargetsUsed,
				minIntensity,
				-1.0f,
				Optional.ofNullable((PeptideModification)null),
				ScoringBreadthType.ENTIRE_RT_WINDOW,
				numberOfExtraDecoyLibrariesSearched,
				quantifyAcrossSamples,
				verifyModificationIons,
				-1.0f,
				SearchParameters.DEFAULT_MIN_NUM_INTEGRATED_RT_POINTS,
				false,
				false, 
				Optional.empty(),
				Optional.empty(),
				true,
				true,
				false,
				false,
				false,
				instrument,
				false
		);
		minPeptideLength=5;
		maxPeptideLength=100;
		this.maxMissedCleavages=maxMissedCleavages;
		this.minCharge=minCharge;
		this.maxCharge=maxCharge;
		numberOfReportedPeaks=1;
		addDecoysToBackgound=false;
		dontRunDecoys=false;
		alpha=1.8f;
		beta=0.4f;
		this.requireVariableMods=requireVariableMods;
	}

	public PecanSearchParameters(
			AminoAcidConstants aaConstants,
			FragmentationType fragType,
			MassTolerance fragmentTolerance,
			MassTolerance precursorTolerance,
			DigestionEnzyme enzyme,
			DataAcquisitionType dataAcquisitionType,
			boolean quantifyAcrossSamples,
			boolean verifyModificationIons,
			boolean requireVariableMods
	) {
		super(
				aaConstants,
				fragType,
				precursorTolerance,
				0.0,
				0.0,
				fragmentTolerance,
				0.0,
				fragmentTolerance,
				enzyme,
				0.01f,
				0.01f,
				false,
				null,
				PercolatorExecutor.DEFAULT_TRAINING_SET_SIZE,
				PercolatorExecutor.DEFAULT_TRAINING_THRESHOLD,
				PercolatorExecutor.DEFAULT_TRAINING_ITERATIONS,
				dataAcquisitionType,
				Runtime.getRuntime().availableProcessors(),
				24f,
				-1f,
				-1f,
				5,
				3,
				-1,
				-1.0f,
				-1.0f,
				Optional.ofNullable((PeptideModification)null),
				ScoringBreadthType.ENTIRE_RT_WINDOW,
				0,
				quantifyAcrossSamples,
				verifyModificationIons,
				-1.0f,
				SearchParameters.DEFAULT_MIN_NUM_INTEGRATED_RT_POINTS,
				false,
				false,
				Optional.empty(),
				Optional.empty(),
				true,
				true,
				false,
				false,
				false,
				InstrumentSpecificSearchParameters.OrbitrapOrbitrap,
				false);
		minPeptideLength=5;
		maxPeptideLength=100;
		maxMissedCleavages=1;
		minCharge=2;
		maxCharge=3;
		numberOfReportedPeaks=1;
		addDecoysToBackgound=false;
		dontRunDecoys=false;
		alpha=1.8f;
		beta=0.4f;
		this.requireVariableMods=requireVariableMods;
	}

	/**
	 * Used only for testing.
	 */
	public PecanSearchParameters(
			AminoAcidConstants aaConstants,
			FragmentationType fragType,
			MassTolerance fragmentTolerance,
			MassTolerance precursorTolerance,
			DigestionEnzyme enzyme,
			boolean quantifyAcrossSamples,
			boolean verifyModificationIons,
			boolean requireVariableMods
	) {
		super(
				aaConstants,
				fragType,
				precursorTolerance,
				0.0,
				0.0,
				fragmentTolerance,
				0.0,
				fragmentTolerance,
				enzyme,
				0.01f,
				0.01f,
				false,
				null,
				PercolatorExecutor.DEFAULT_TRAINING_SET_SIZE,
				PercolatorExecutor.DEFAULT_TRAINING_THRESHOLD,
				PercolatorExecutor.DEFAULT_TRAINING_ITERATIONS,
				DataAcquisitionType.DIA,
				Runtime.getRuntime().availableProcessors(),
				24f,
				-1f,
				-1f,
				5,
				3,
				-1,
				-1.0f,
				-1.0f,
				Optional.ofNullable((PeptideModification)null),
				ScoringBreadthType.ENTIRE_RT_WINDOW,
				0,
				quantifyAcrossSamples,
				verifyModificationIons,
				-1.0f,
				SearchParameters.DEFAULT_MIN_NUM_INTEGRATED_RT_POINTS,
				false,
				false, 
				Optional.empty(),
				Optional.empty(),
				true,
				true,
				false,
				false,
				false,
				InstrumentSpecificSearchParameters.OrbitrapOrbitrap,
				false
		);
		minPeptideLength=5;
		maxPeptideLength=100;
		maxMissedCleavages=1;
		minCharge=2;
		maxCharge=3;
		numberOfReportedPeaks=1;
		addDecoysToBackgound=false;
		dontRunDecoys=false;
		alpha=1.8f;
		beta=0.4f;
		this.requireVariableMods=requireVariableMods;
	}

	/**
	 * Used only for testing.
	 */
	public PecanSearchParameters(
			AminoAcidConstants aaConstants,
			FragmentationType fragType,
			MassTolerance fragmentTolerance,
			MassTolerance precursorTolerance,
			DigestionEnzyme enzyme,
			int maxMissedCleavages,
			boolean quantifyAcrossSamples,
			boolean verifyModificationIons,
			boolean requireVariableMods
	) {
		super(
				aaConstants,
				fragType,
				precursorTolerance,
				0.0,
				0.0,
				fragmentTolerance,
				0.0,
				fragmentTolerance,
				enzyme,
				0.01f,
				0.01f,
				false,
				null,
				PercolatorExecutor.DEFAULT_TRAINING_SET_SIZE,
				PercolatorExecutor.DEFAULT_TRAINING_THRESHOLD,
				PercolatorExecutor.DEFAULT_TRAINING_ITERATIONS,
				DataAcquisitionType.DIA,
				Runtime.getRuntime().availableProcessors(),
				24f,
				-1f,
				-1f,
				5,
				3,
				-1,
				-1.0f,
				-1.0f,
				Optional.ofNullable((PeptideModification)null),
				ScoringBreadthType.ENTIRE_RT_WINDOW,
				0,
				quantifyAcrossSamples,
				verifyModificationIons,
				-1.0f,
				SearchParameters.DEFAULT_MIN_NUM_INTEGRATED_RT_POINTS,
				false,
				false, 
				Optional.empty(),
				Optional.empty(),
				true,
				true,
				false,
				false,
				false,
				InstrumentSpecificSearchParameters.OrbitrapOrbitrap,
				false);
		this.maxMissedCleavages=maxMissedCleavages;
		minPeptideLength=5;
		maxPeptideLength=100;
		minCharge=2;
		maxCharge=3;
		numberOfReportedPeaks=1;
		addDecoysToBackgound=false;
		dontRunDecoys=false;
		alpha=1.8f;
		beta=0.4f;
		this.requireVariableMods=requireVariableMods;
	}

	public int getMaxMissedCleavages() {
		return maxMissedCleavages;
	}

	public int getMaxPeptideLength() {
		return maxPeptideLength;
	}

	public int getMinPeptideLength() {
		return minPeptideLength;
	}

	public byte getMaxCharge() {
		return maxCharge;
	}

	public byte getMinCharge() {
		return minCharge;
	}

	public int getNumberOfReportedPeaks() {
		return numberOfReportedPeaks;
	}

	public boolean isAddDecoysToBackgound() {
		return addDecoysToBackgound;
	}
	
	public boolean isDontRunDecoys() {
		return dontRunDecoys;
	}

	public boolean isHandleNTermMethionineInDigestion() {
		return PecanSearchParameters.DEFAULT_HANDLE_N_TERM_METHIONINE;
	}

	public float getAlpha() {
		return alpha;
	}
	
	public float getBeta() {
		return beta;
	}
	public boolean isRequireVariableMods() {
		return requireVariableMods;
	}
}
