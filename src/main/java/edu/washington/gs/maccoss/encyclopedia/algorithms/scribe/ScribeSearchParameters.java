package edu.washington.gs.maccoss.encyclopedia.algorithms.scribe;

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

import edu.washington.gs.maccoss.encyclopedia.Encyclopedia;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorVersion;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ScoringBreadthType;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ThesaurusSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

public class ScribeSearchParameters extends SearchParameters {
	@Override
	public void savePreferences(File libraryFile, File fastaFile) throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("Scribe");
		HashMap<String, String> map=toParameterMap();
		if (libraryFile!=null) map.put(Encyclopedia.TARGET_LIBRARY_TAG, libraryFile.getAbsolutePath());
		if (fastaFile!=null) map.put(Encyclopedia.BACKGROUND_FASTA_TAG, fastaFile.getAbsolutePath());
		for (Entry<String, String> entry : map.entrySet()) {
			//System.out.println("Writing Scribe preference "+entry.getKey()+" = "+entry.getValue());
			prefs.put(entry.getKey(), entry.getValue());
		}
		prefs.flush();
	}
	
	public static HashMap<String, String> readPreferences() throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("Scribe");
		HashMap<String, String> map=new HashMap<String, String>();
		for (String key : prefs.keys()) {
			String value=prefs.get(key, "");
			//System.out.println("Reading Scribe preference "+key+" = "+value);
			map.put(key, value);
		}
		return map;
	}

	@Override
	public HashMap<String, String> toParameterMap() {
		HashMap<String, String> map=super.toParameterMap();

		return map;
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
	
	public static HashMap<String, String> getDefaultParameters() {
		HashMap<String, String> defaults=SearchParameterParser.getDefaultParameters();
		return defaults;
	}
	
	public static ScribeSearchParameters parseParameters(HashMap<String, String> map) {
		SearchParameters params=SearchParameterParser.parseParameters(map);
		return convertFromEncyclopeDIA(params);
	}
	
	public static SearchParameters getDefaultParametersObject() {
		return parseParameters(getDefaultParameters());
	}
	
	public static ScribeSearchParameters readFromXML(Document doc, Element rootElement) {
		if (!rootElement.getTagName().equals(ScribeSearchParameters.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+ScribeSearchParameters.class.getSimpleName()+"]");
		}
		
		HashMap<String, String> paramMap=getDefaultParametersObject().toParameterMap();
		
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
		return parseParameters(paramMap);
	}


	public ScribeSearchParameters(
			AminoAcidConstants aaConstants,
			FragmentationType fragType,
			MassTolerance precursorTolerance,
			double precursorOffsetPPM,
			double precursorIsolationMargin,
			MassTolerance fragmentTolerance,
			double fragmentOffsetPPM,
			MassTolerance libraryFragmentTolerance,
			DigestionEnzyme enzyme,
			float percolatorThreshold,
			float percolatorProteinThreshold,
			PercolatorVersion percolatorVersionNumber,
			int percolatorTrainingSetSize,
			float percolatorTrainingSetThreshold,
			int percolatorTrainingIterations,
			DataAcquisitionType dataAcquisitionType,
			int numberOfThreadsUsed,
			float expectedPeakWidth,
			float targetWindowCenter,
			float precursorWindowSize,
			int numberOfQuantitativePeaks,
			int minNumOfQuantitativePeaks,
			int topNTargetsUsed,
			float minIntensity,
			float minIntensityNumIons, 
			Optional<PeptideModification> modification,
			ScoringBreadthType searchType,
			float getNumberOfExtraDecoyLibrariesSearched,
			boolean quantifyAcrossSamples,
			boolean verifyModificationIons,
			boolean filterPeaklists,
			boolean doNotUseGlobalFDR,
			Optional<File> precursorIsolationRangeFile, 
			Optional<File> percolatorModelFile,
			boolean normalizeByTIC,
			boolean subtractBackground,
			boolean maskBadIntegrations,
			boolean integratePrecursors,
			boolean enableAdvancedOptions
	) {
		super(
				aaConstants,
				fragType,
				precursorTolerance,
				precursorOffsetPPM,
				precursorIsolationMargin,
				fragmentTolerance,
				fragmentOffsetPPM,
				libraryFragmentTolerance,
				enzyme,
				percolatorThreshold,
				percolatorProteinThreshold,
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
				modification,
				searchType,
				getNumberOfExtraDecoyLibrariesSearched,
				quantifyAcrossSamples,
				verifyModificationIons,
				-1.0f,
				filterPeaklists,
				doNotUseGlobalFDR,
				precursorIsolationRangeFile,
				percolatorModelFile,
				normalizeByTIC,
				subtractBackground,
				maskBadIntegrations,
				integratePrecursors,
				enableAdvancedOptions
		);
	}

	public static ScribeSearchParameters convertFromEncyclopeDIA(SearchParameters params) {
		return new ScribeSearchParameters(
				params.getAAConstants(),
				params.getFragType(),
				params.getPrecursorTolerance(),
				params.getPrecursorOffsetPPM(),
				params.getPrecursorIsolationMargin(),
				params.getFragmentTolerance(),
				params.getFragmentOffsetPPM(),
				params.getLibraryFragmentTolerance(),
				params.getEnzyme(),
				params.getPercolatorThreshold(),
				params.getPercolatorProteinThreshold(),
				params.getPercolatorVersionNumber(),
				params.getPercolatorTrainingSetSize(),
				params.getPercolatorTrainingSetThreshold(),
				params.getPercolatorTrainingIterations(),
				params.getDataAcquisitionType(),
				params.getNumberOfThreadsUsed(),
				params.getExpectedPeakWidth(),
				params.getTargetWindowCenter(),
				params.getPrecursorWindowSize(),
				params.getNumberOfQuantitativePeaks(),
				params.getMinNumOfQuantitativePeaks(),
				params.getTopNTargetsUsed(),
				params.getMinIntensity(),
				params.getMinIntensityNumIons(),
				params.getLocalizingModification(),
				params.getScoringBreadthType(),
				params.getNumberOfExtraDecoyLibrariesSearched(),
				params.isQuantifySameFragmentsAcrossSamples(),
				params.isVerifyModificationIons(),
				params.isFilterPeaklists(),
				params.isDoNotUseGlobalFDR(),
				params.getPrecursorIsolationRangeFile(),
				params.getPercolatorModelFile(),
				params.isNormalizeByTIC(),
				params.isSubtractBackground(),
				params.isMaskBadIntegrations(),
				params.isIntegratePrecursors(),
				params.isEnableAdvancedOptions()
		);
	}
}
