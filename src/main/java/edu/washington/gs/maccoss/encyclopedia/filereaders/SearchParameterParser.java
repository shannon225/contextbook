package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.HashMap;
import java.util.NoSuchElementException;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Iterables;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorVersion;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ScoringBreadthType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.parameters.InstrumentSpecificSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassErrorUnitType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class SearchParameterParser {

	public static HashMap<String,String> getDefaultParameters() {
		HashMap<String, String> map=new HashMap<String, String>();
		map.put("-fixed", "C=57.0214635");
		map.put("-frag", "CID");
		map.put("-ptol", "10");
		map.put("-ftol", "10");
		map.put("-lftol", "10");
		map.put("-ptolunits", "ppm");
		map.put("-ftolunits", "ppm");
		map.put("-lftolunits", "ppm");
		map.put("-poffset", "0");
		map.put("-foffset", "0");
		map.put("-precursorIsolationMargin", "0");
		map.put("-precursorWindowSize", "-1");
		map.put("-enzyme", "trypsin");
		map.put("-usePercolator", "true");
		map.put("-percolatorThreshold", "0.01");
		map.put("-percolatorProteinThreshold", "0.01");
		map.put("-percolatorVersion", PercolatorExecutor.DEFAULT_VERSION_NUMBER.toString());
		map.put(SearchParameters.OPT_PERC_TRAINING_SIZE, Integer.toString(PercolatorExecutor.DEFAULT_TRAINING_SET_SIZE));
		map.put(SearchParameters.OPT_PERC_TRAINING_THRESH, Float.toString(PercolatorExecutor.DEFAULT_TRAINING_THRESHOLD));
		map.put("-expectedPeakWidth", "25");
		map.put("-maxWindowWidth", "-1");
		map.put("-acquisition", DataAcquisitionType.toString(DataAcquisitionType.DIA));
		map.put("-localizationModification", PeptideModification.NO_MODIFICATION_NAME);
		map.put("-scoringBreadthType", ScoringBreadthType.RECALIBRATED_PEAK_WIDTH.toShortname());
		map.put("-numberOfExtraDecoyLibrariesSearched", "0.0");
		map.put(SearchParameters.NUMBER_OF_QUANTITATIVE_PEAKS, "5");
		map.put("-minNumOfQuantitativePeaks", "3");
		map.put("-topNTargetsUsed", ""+SearchParameters.DEFAULT_MAX_NUM_FRAGMENT_IONS);
		map.put("-verifyModificationIons", "true");
		map.put("-minIntensity", "-1.0");
		map.put("-minIntensityNumIons", "-1.0");
		map.put("-rtWindowInMin", "-1.0");
        map.put("-filterPeaklists", "false");
        map.put("-numberOfThreadsUsed", Integer.toString(Runtime.getRuntime().availableProcessors()));
		map.put("-normalizeByTIC", "true");
		map.put(SearchParameters.SUBTRACT_BACKGROUND, "true");
		map.put(SearchParameters.SMOOTH_INTEGRATIONS, "true");
		map.put(SearchParameters.MASK_BAD_INTEGRATIONS, "false");
		map.put(SearchParameters.INTEGRATE_PRECURSORS, "false");
		map.put(SearchParameters.ADJUST_INFERRED_RT_BOUNDARIES, "false");
		map.put(SearchParameters.MIN_NUM_INTEGRATED_RT_POINTS, Integer.toString(SearchParameters.DEFAULT_MIN_NUM_INTEGRATED_RT_POINTS));
		map.put("-skipLibraryRetentionTime", "false");
		return map;
	}
	
	/**
	 * parameters that can affect file exports
	 * @return
	 */
	public static HashMap<String,String> getExportParameters() {
		HashMap<String, String> map=new HashMap<String, String>();
		map.put("-fixed", "C=57.0214635");
		map.put("-ftol", "10");
		map.put("-ftolunits", "ppm");
		map.put("-foffset", "0");
		map.put("-percolatorThreshold", "0.01");
		map.put("-percolatorProteinThreshold", "0.01");
		map.put("-percolatorLocation", "internal");
		map.put("-localizationModification", PeptideModification.NO_MODIFICATION_NAME);
		map.put("-numberOfExtraDecoyLibrariesSearched", "0.0");
		map.put(SearchParameters.NUMBER_OF_QUANTITATIVE_PEAKS, "5");
		map.put("-minNumOfQuantitativePeaks", "3");
		map.put("-normalizeByTIC", "true");
		map.put(SearchParameters.SUBTRACT_BACKGROUND, "true");
		map.put(SearchParameters.SMOOTH_INTEGRATIONS, "true");
		map.put(SearchParameters.MASK_BAD_INTEGRATIONS, "false");
		map.put(SearchParameters.INTEGRATE_PRECURSORS, "false");
		return map;
	}
	
	public static SearchParameters getDefaultParametersObject(InstrumentSpecificSearchParameters instrument) {
		HashMap<String, String> defaultParameters = instrument.overwriteParameters(getDefaultParameters());
		return parseParameters(defaultParameters);
	}
	
	public static SearchParameters getDefaultParametersObject() {
		return parseParameters(getDefaultParameters());
	}
	
	public static SearchParameters parseParameters(File defaultParameters, HashMap<String, String> parameters) {
		HashMap<String, String> map=ParsingUtils.readFile(defaultParameters);
		map.putAll(parameters);
		return parseParameters(map);
	}
	
	public static SearchParameters parseParameters(HashMap<String, String> parameters) {
		final AminoAcidConstants aaConstants=new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());
		final FragmentationType fragType;
		
		final MassErrorUnitType precursorToleranceType;
		final MassErrorUnitType fragmentToleranceType;
		final MassErrorUnitType libraryFragmentToleranceType;
		final MassTolerance precursorTolerance;
		final MassTolerance fragmentTolerance;
		final MassTolerance libraryFragmentTolerance;
		final double precursorOffsetPPM;
		final double fragmentOffsetPPM;
		final double precursorIsolationMargin;
		final DigestionEnzyme enzyme;
		final float percolatorThreshold;
		final float percolatorProteinThreshold;
		final boolean usePercolator;
		final PercolatorVersion percolatorVersionNumber;
		final int percolatorTrainingSetSize;
		final float percolatorTrainingSetThreshold;
		final int percolatorTrainingIterations;
		final DataAcquisitionType dataAcquisitionType;
		final int numberOfThreadsUsed;
		final float targetWindowCenter;
		final float expectedPeakWidth;
		final float precursorWindowSize;
		final float maxWindowWidth;
		final int numberOfQuantitativePeaks;
		final int minNumOfQuantitativePeaks;
		final int topNTargetsUsed;
		final float minIntensity;
		final float minIntensityNumIons;
		final float numberOfExtraDecoyLibrariesSearched;
		final Optional<PeptideModification> localizationModification;
		final ScoringBreadthType breadthType;
		final boolean quantifyAcrossSamples;
		final boolean verifyModificationIons;
		final float rtWindowInMin;
		final int minNumIntegratedRTPoints;
        final boolean filterPeaklists;
        final boolean doNotUseGlobalFDR;
        final boolean normalizeByTIC;
        final boolean subtractBackground;
        final boolean smoothIntegrations;
        final boolean maskBadIntegrations;
        final boolean integratePrecursors;
        final boolean adjustInferredRTBoundaries;
        final boolean enableAdvancedOptions;
        final boolean skipLibraryRetentionTime;
        final Optional<File> percolatorModelFile;
        final Optional<File> precursorIsolationRangeFile;
		
		String value=parameters.get("-frag");
		if (value==null) {
			fragType=FragmentationType.CID;
		} else {
			fragType=FragmentationType.getFragmentationType(value);
		}
		if (fragType==null) {
			throw new EncyclopediaException("Error parsing fragmentation type from ["+value+"]");
		}
		
		value=parameters.get("-acquisition");
		if (value==null) {
			dataAcquisitionType=DataAcquisitionType.DIA;
		} else {
			dataAcquisitionType=DataAcquisitionType.getAcquisitionType(value);
		}
		if (dataAcquisitionType==null) {
			throw new EncyclopediaException("Error parsing acquisition type from ["+value+"]");
		}
		
		value=parameters.get("-ptolunits");
		if (value==null) {
			precursorToleranceType=MassErrorUnitType.PPM;
		} else {
			precursorToleranceType=MassErrorUnitType.getUnitType(value);
		}
		if (precursorToleranceType==null) {
			throw new EncyclopediaException("Error parsing precursor mass error unit type from ["+value+"]");
		}
		
		value=parameters.get("-ftolunits");
		if (value==null) {
			fragmentToleranceType=MassErrorUnitType.PPM;
		} else {
			fragmentToleranceType=MassErrorUnitType.getUnitType(value);
		}
		if (fragmentToleranceType==null) {
			throw new EncyclopediaException("Error parsing fragment mass error unit type from ["+value+"]");
		}
		
		value=parameters.get("-lftolunits");
		if (value==null) {
			libraryFragmentToleranceType=MassErrorUnitType.PPM;
		} else {
			libraryFragmentToleranceType=MassErrorUnitType.getUnitType(value);
		}
		if (libraryFragmentToleranceType==null) {
			throw new EncyclopediaException("Error parsing library mass error unit type from ["+value+"]");
		}
		
		value=parameters.get("-ptol");
		if (value==null) {
			precursorTolerance=new MassTolerance(10, MassErrorUnitType.PPM);
		} else {
			try {
				precursorTolerance=new MassTolerance(Double.parseDouble(value), precursorToleranceType);
			} catch (NumberFormatException nfe) {
				throw new EncyclopediaException("Error parsing precursor tolerance from ["+value+"]", nfe);
			}
		}
		
		value=parameters.get("-ftol");
		if (value==null) {
			fragmentTolerance=new MassTolerance(10, MassErrorUnitType.PPM);
		} else {
			try {
				fragmentTolerance=new MassTolerance(Double.parseDouble(value), fragmentToleranceType);
			} catch (NumberFormatException nfe) {
				throw new EncyclopediaException("Error parsing fragment tolerance from ["+value+"]", nfe);
			}
		}
		
		value=parameters.get("-lftol");
		if (value==null) {
			libraryFragmentTolerance=new MassTolerance(10, MassErrorUnitType.PPM);
		} else {
			try {
				libraryFragmentTolerance=new MassTolerance(Double.parseDouble(value), libraryFragmentToleranceType);
			} catch (NumberFormatException nfe) {
				throw new EncyclopediaException("Error parsing library fragment tolerance from ["+value+"]", nfe);
			}
		}
		
		value=parameters.get("-poffset");
		if (value==null) {
			precursorOffsetPPM=0.0;
		} else {
			try {
				precursorOffsetPPM=Double.parseDouble(value);
			} catch (NumberFormatException nfe) {
				throw new EncyclopediaException("Error parsing precursor tolerance from ["+value+"]", nfe);
			}
		}
		
		value=parameters.get("-foffset");
		if (value==null) {
			fragmentOffsetPPM=0.0;
		} else {
			try {
				fragmentOffsetPPM=Double.parseDouble(value);
			} catch (NumberFormatException nfe) {
				throw new EncyclopediaException("Error parsing fragment tolerance from ["+value+"]", nfe);
			}
		}
		
		value=parameters.get("-precursorIsolationMargin");
		if (value==null) {
			precursorIsolationMargin=0.0;
		} else {
			try {
				precursorIsolationMargin=Double.parseDouble(value);
			} catch (NumberFormatException nfe) {
				throw new EncyclopediaException("Error parsing precursor isolation margin from ["+value+"]", nfe);
			}
		}
		
		value=parameters.get("-enzyme");
		if (value==null) {
			enzyme=DigestionEnzyme.getEnzyme("trypsin");
		} else {
			enzyme=DigestionEnzyme.getEnzyme(value);
		}
		value=parameters.get("-percolatorModelFile");
		if (value==null) {
			percolatorModelFile=Optional.empty();
		} else {
			File f=new File(value);
			if (f.exists()&&f.canRead()) {
				percolatorModelFile=Optional.of(f);
			} else {
				percolatorModelFile=Optional.empty();
			}
		}
		value=parameters.get("-precursorIsolationRangeFile");
		if (value==null) {
			precursorIsolationRangeFile=Optional.empty();
		} else {
			File f=new File(value);
			if (f.exists()&&f.canRead()) {
				precursorIsolationRangeFile=Optional.of(f);
			} else {
				precursorIsolationRangeFile=Optional.empty();
			}
		}
		
		usePercolator=ParsingUtils.getBoolean("-usePercolator", parameters, true);
		percolatorThreshold=ParsingUtils.getFloat("-percolatorThreshold", parameters, 0.01f);
		percolatorProteinThreshold=ParsingUtils.getFloat("-percolatorProteinThreshold", parameters, 0.01f);
		percolatorVersionNumber=PercolatorVersion.getVersion(parameters.get("-percolatorVersion"));
		percolatorTrainingSetSize = ParsingUtils.getInteger(SearchParameters.OPT_PERC_TRAINING_SIZE, parameters, PercolatorExecutor.DEFAULT_TRAINING_SET_SIZE);
		percolatorTrainingSetThreshold = ParsingUtils.getFloat(SearchParameters.OPT_PERC_TRAINING_THRESH, parameters, PercolatorExecutor.DEFAULT_TRAINING_THRESHOLD);
		percolatorTrainingIterations= ParsingUtils.getInteger("-percolatorTrainingIterations", parameters, PercolatorExecutor.DEFAULT_TRAINING_ITERATIONS);

		numberOfThreadsUsed=ParsingUtils.getInteger("-numberOfThreadsUsed", parameters, Runtime.getRuntime().availableProcessors());
		targetWindowCenter=ParsingUtils.getFloat("-targetWindowCenter", parameters, -1f);
		precursorWindowSize=ParsingUtils.getFloat("-precursorWindowSize", parameters, -1f);
		maxWindowWidth=ParsingUtils.getFloat("-maxWindowWidth", parameters, -1f);
		expectedPeakWidth=ParsingUtils.getFloat("-expectedPeakWidth", parameters, 25f);
		numberOfQuantitativePeaks=ParsingUtils.getInteger(SearchParameters.NUMBER_OF_QUANTITATIVE_PEAKS, parameters, 5);
		minNumOfQuantitativePeaks=ParsingUtils.getInteger("-minNumOfQuantitativePeaks", parameters, 3);
		topNTargetsUsed=ParsingUtils.getInteger("-topNTargetsUsed", parameters, SearchParameters.DEFAULT_MAX_NUM_FRAGMENT_IONS);
		
		minIntensity=ParsingUtils.getFloat("-minIntensity", parameters, -1.0f);
		minIntensityNumIons=ParsingUtils.getFloat("-minIntensityNumIons", parameters, -1.0f);
		rtWindowInMin=ParsingUtils.getFloat("-rtWindowInMin", parameters, -1f);
		minNumIntegratedRTPoints=ParsingUtils.getInteger(SearchParameters.MIN_NUM_INTEGRATED_RT_POINTS, parameters, SearchParameters.DEFAULT_MIN_NUM_INTEGRATED_RT_POINTS);
		
		value=parameters.get("-localizationModification");
		if (value != null) {
			final String localizationModificationName = value;
			Set<PeptideModification> peptideModifications =
					aaConstants.getLocalizationModifications()
							.stream()
							.filter(mod -> localizationModificationName.equalsIgnoreCase(mod.getShortname()))
							.collect(Collectors.toSet());

			PeptideModification mod;
			try {
				mod = Iterables.getOnlyElement(peptideModifications);
			} catch (NoSuchElementException noElement) {
				// Preserves previous behavior where a mod 'unknown' to the system will be treated as not specifying a localization mod.
				// We think we should throw in this case since this silent error is misleading.
				mod = null;
			} catch (IllegalStateException multipleElements) {
				throw new IllegalStateException("Multiple modifications correspond to " + localizationModificationName);
			}

			localizationModification=Optional.ofNullable(mod);

		} else {
			localizationModification=Optional.empty();
		}
		
		value=parameters.get("-scoringBreadthType");
		if (value!=null) {
			ScoringBreadthType type;
			try {
				type=ScoringBreadthType.getType(value);
			} catch (Exception e) {
				Logger.errorLine("Falling back to scoring breadth type: "+ScoringBreadthType.RECALIBRATED_PEAK_WIDTH.toShortname());
				type=ScoringBreadthType.RECALIBRATED_PEAK_WIDTH;
			}
			breadthType=type;
		} else {
			breadthType=ScoringBreadthType.RECALIBRATED_PEAK_WIDTH;
		}
		
		float tempNumberOfExtraDecoyLibrariesSearched=ParsingUtils.getFloat("-numberOfExtraDecoyLibrariesSearched", parameters, 0.0f);
		if (tempNumberOfExtraDecoyLibrariesSearched<0.0f) {
			Logger.errorLine("-numberOfExtraDecoyLibrariesSearched cannot be less than 0%! Using 0% extra decoys.");
			numberOfExtraDecoyLibrariesSearched=0.0f;
		} else {
			numberOfExtraDecoyLibrariesSearched=tempNumberOfExtraDecoyLibrariesSearched;
		}
		quantifyAcrossSamples=ParsingUtils.getBoolean("-quantifyAcrossSamples", parameters, true); // edited to match GUI, see PeptideQuantExtractorTask line 252
		verifyModificationIons=ParsingUtils.getBoolean("-verifyModificationIons", parameters, true);
        filterPeaklists=ParsingUtils.getBoolean("-filterPeaklists", parameters, false);
        doNotUseGlobalFDR=ParsingUtils.getBoolean("-doNotUseGlobalFDR", parameters, false);
		normalizeByTIC = ParsingUtils.getBoolean("-normalizeByTIC", parameters, true);
        subtractBackground=ParsingUtils.getBoolean(SearchParameters.SUBTRACT_BACKGROUND, parameters, true);
        smoothIntegrations=ParsingUtils.getBoolean(SearchParameters.SMOOTH_INTEGRATIONS, parameters, true);
        maskBadIntegrations=ParsingUtils.getBoolean(SearchParameters.MASK_BAD_INTEGRATIONS, parameters, false);
        integratePrecursors=ParsingUtils.getBoolean(SearchParameters.INTEGRATE_PRECURSORS, parameters, false);
        adjustInferredRTBoundaries=ParsingUtils.getBoolean(SearchParameters.ADJUST_INFERRED_RT_BOUNDARIES, parameters, false);
        enableAdvancedOptions=ParsingUtils.getBoolean(SearchParameters.ENABLE_ADVANCED_OPTIONS, parameters, false);
        skipLibraryRetentionTime=ParsingUtils.getBoolean("-skipLibraryRetentionTime", parameters, false);
		 
        String instrumentValue=parameters.get(SearchParameters.INSTRUMENT);
        InstrumentSpecificSearchParameters instrument=instrumentValue==null?InstrumentSpecificSearchParameters.OrbitrapOrbitrap:InstrumentSpecificSearchParameters.fromString(instrumentValue);

		return new SearchParameters(
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
				maxWindowWidth,
				numberOfQuantitativePeaks,
				minNumOfQuantitativePeaks,
				topNTargetsUsed,
				minIntensity,
				minIntensityNumIons,
				localizationModification,
				breadthType,
				numberOfExtraDecoyLibrariesSearched,
				quantifyAcrossSamples,
				verifyModificationIons,
				rtWindowInMin,
				minNumIntegratedRTPoints,
				filterPeaklists,
				doNotUseGlobalFDR,
				precursorIsolationRangeFile,
				percolatorModelFile,
				smoothIntegrations,
				normalizeByTIC,
				subtractBackground,
				maskBadIntegrations,
				adjustInferredRTBoundaries,
				skipLibraryRetentionTime,
				integratePrecursors,
				instrument,
				enableAdvancedOptions
		);
	}
}
