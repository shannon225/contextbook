package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.HashMap;
import java.util.Optional;
import java.util.StringTokenizer;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorVersion;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassErrorUnitType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class PecanParameterParser {
	public static HashMap<String,String> getDefaultParameters() {
		HashMap<String, String> map=new HashMap<String, String>();
		map.put("-fixed", "C=57.0214635");
		map.put("-variable", "-");
		map.put("-frag", "YONLY");
		map.put("-ptol", "10");
		map.put("-ftol", "10");
		map.put("-ptolunits", "ppm");
		map.put("-ftolunits", "ppm");
		map.put("-poffset", "0");
		map.put("-foffset", "0");
		map.put("-precursorIsolationMargin", "0");
		map.put("-enzyme", "trypsin");
		map.put("-minLength", "5");
		map.put("-maxLength", "100");
		map.put("-maxMissedCleavage", "1");
		map.put("-minCharge", "2");
		map.put("-maxCharge", "3");
		map.put("-expectedPeakWidth", "25");
		map.put("-numberOfReportedPeaks", "1");
		map.put("-addDecoysToBackground", "false");
		map.put("-dontRunDecoys", "false");
		map.put("-percolatorThreshold", "0.01");
		map.put("-percolatorProteinThreshold", "0.01");
		map.put("-alpha", "1.8");
		map.put("-beta", "0.4");
		map.put("-percolatorVersion", PercolatorExecutor.DEFAULT_VERSION_NUMBER.toString());
		map.put("-acquisition", DataAcquisitionType.toString(DataAcquisitionType.DIA));
		map.put("-precursorWindowSize", "-1");
		map.put("-numberOfThreadsUsed", Integer.toString(Runtime.getRuntime().availableProcessors()));
		map.put(SearchParameters.NUMBER_OF_QUANTITATIVE_PEAKS, "5");
		map.put("-minNumOfQuantitativePeaks", "3");
		map.put("-minIntensity", "-1.0");
		map.put("-requireVariableMods", "false");
        map.put("-filterPeaklists", "false");
        map.put("-normalizeByTIC", "true");
		return map;
	}
	
	public static PecanSearchParameters getDefaultParametersObject() {
		return parseParameters(getDefaultParameters());
	}
	
	public static PecanSearchParameters parseParameters(File defaultParameters, HashMap<String, String> parameters) {
		HashMap<String, String> map=ParsingUtils.readFile(defaultParameters);
		map.putAll(parameters);
		return parseParameters(map);
	}
	
	public static PecanSearchParameters parseParameters(HashMap<String, String> parameters) {
		final AminoAcidConstants aaConstants;
		final FragmentationType fragType;
		final MassErrorUnitType precursorToleranceType;
		final MassErrorUnitType fragmentToleranceType;
		final MassTolerance precursorTolerance;
		final MassTolerance fragmentTolerance;
		final double precursorOffsetPPM;
		final double fragmentOffsetPPM;
		final double precursorIsolationMargin;
		final DigestionEnzyme enzyme;
		final int minPeptideLength;
		final int maxPeptideLength;
		final int maxMissedCleavages;
		final byte minCharge;
		final byte maxCharge;
		final float expectedPeakWidth;
		final int numberOfReportedPeaks;
		final boolean addDecoysToBackgound;
		final boolean dontRunDecoys;
		final float percolatorThreshold;
		final float percolatorProteinThreshold;
		final PercolatorVersion percolatorVersionNumber;
		final int percolatorTrainingSetSize;
		final float percolatorTrainingSetThreshold;
		final int percolatorTrainingIterations;
		final float alpha;
		final float beta;
		final DataAcquisitionType dataAcquisitionType;
		final int numberOfThreadsUsed;
		final float targetWindowCenter;
		final float precursorWindowSize;
		final int numberOfQuantitativePeaks;
		final int minNumOfQuantitativePeaks;
		final int topNTargetsUsed;
		final float minIntensity;
		final boolean quantifyAcrossSamples;
		final boolean requireVariableMods;
        final boolean filterPeaklists;
        final boolean doNotUseGlobalFDR;
        final Optional<File> percolatorModelFile;
        final Optional<File> precursorIsolationRangeFile;
        final boolean normalizeByTIC;
        final boolean subtractBackground;
        final boolean enableAdvancedOptions;

		ModificationMassMap variableMods=new ModificationMassMap(parameters.get("-variable"));

		TCharDoubleHashMap fixedMods=new TCharDoubleHashMap();
		String value=parameters.get("-fixed");
		if (value!=null) {
			try {
				StringTokenizer st=new StringTokenizer(value, ",");
				while (st.hasMoreTokens()) {
					String token=st.nextToken();
					char aa=token.charAt(0);
					double mass=Double.parseDouble(token.substring(2)); // +1 for '=' (could actually be any deliminator)
					fixedMods.put(aa, mass);
				}
			} catch (Exception e) {
				throw new EncyclopediaException("Error parsing fixed modifications from ["+value+"]", e);
			}
		} else {
			fixedMods.put('C', 57.0214635);
		}
		aaConstants=new AminoAcidConstants(fixedMods, variableMods);

		value=parameters.get("-frag");
		if (value==null) {
			fragType=FragmentationType.HCD;
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

		minPeptideLength=ParsingUtils.getInteger("-minLength", parameters, 5);
		maxPeptideLength=ParsingUtils.getInteger("-maxLength", parameters, 100);
		maxMissedCleavages=ParsingUtils.getInteger("-maxMissedCleavage", parameters, 1);
		minCharge=ParsingUtils.getByte("-minCharge", parameters, (byte)2);
		maxCharge=ParsingUtils.getByte("-maxCharge", parameters, (byte)3);
		expectedPeakWidth=ParsingUtils.getFloat("-expectedPeakWidth", parameters, 25.0f);
		numberOfReportedPeaks=ParsingUtils.getInteger("-numberOfReportedPeaks", parameters, 1);
		addDecoysToBackgound=ParsingUtils.getBoolean("-addDecoysToBackground", parameters, false);
		dontRunDecoys=ParsingUtils.getBoolean("-dontRunDecoys", parameters, false);

		percolatorThreshold=ParsingUtils.getFloat("-percolatorThreshold", parameters, 0.01f);
		percolatorProteinThreshold=ParsingUtils.getFloat("-percolatorProteinThreshold", parameters, 0.01f);
		percolatorVersionNumber=PercolatorVersion.getVersion(parameters.get("-percolatorVersion"));
		percolatorTrainingSetSize = ParsingUtils.getInteger(SearchParameters.OPT_PERC_TRAINING_SIZE, parameters, PercolatorExecutor.DEFAULT_TRAINING_SET_SIZE);
		percolatorTrainingSetThreshold = ParsingUtils.getFloat(SearchParameters.OPT_PERC_TRAINING_THRESH, parameters, PercolatorExecutor.DEFAULT_TRAINING_THRESHOLD);
		percolatorTrainingIterations= ParsingUtils.getInteger("-percolatorTrainingIterations", parameters, PercolatorExecutor.DEFAULT_TRAINING_ITERATIONS);

		alpha=ParsingUtils.getFloat("-alpha", parameters, 1.8f);
		beta=ParsingUtils.getFloat("-beta", parameters, 0.4f);
		numberOfThreadsUsed=ParsingUtils.getInteger("-numberOfThreadsUsed", parameters, Runtime.getRuntime().availableProcessors());
		targetWindowCenter=ParsingUtils.getFloat("-targetWindowCenter", parameters, -1f);
		precursorWindowSize=ParsingUtils.getFloat("-precursorWindowSize", parameters, -1f);
		numberOfQuantitativePeaks=ParsingUtils.getInteger(SearchParameters.NUMBER_OF_QUANTITATIVE_PEAKS, parameters, 5);
		minNumOfQuantitativePeaks=ParsingUtils.getInteger("-minNumOfQuantitativePeaks", parameters, 3);
		topNTargetsUsed=ParsingUtils.getInteger("-topNTargetsUsed", parameters, -1);
		
		minIntensity=ParsingUtils.getFloat("-minIntensity", parameters, -1.0f);
		quantifyAcrossSamples=ParsingUtils.getBoolean("-quantifyAcrossSamples", parameters, false);
		requireVariableMods=ParsingUtils.getBoolean("-requireVariableMods", parameters, false);
		filterPeaklists = ParsingUtils.getBoolean("-filterPeaklists", parameters, false);
		doNotUseGlobalFDR = ParsingUtils.getBoolean("-doNotUseGlobalFDR", parameters, false);
		subtractBackground = ParsingUtils.getBoolean(SearchParameters.SUBTRACT_BACKGROUND, parameters, true);
		normalizeByTIC = ParsingUtils.getBoolean("-normalizeByTIC", parameters, true);
		enableAdvancedOptions = ParsingUtils.getBoolean(SearchParameters.ENABLE_ADVANCED_OPTIONS, parameters, false);

		return new PecanSearchParameters(
				aaConstants,
				fragType,
				precursorTolerance,
				precursorOffsetPPM,
				precursorIsolationMargin,
				fragmentTolerance,
				fragmentOffsetPPM,
				enzyme,
				expectedPeakWidth,
				minPeptideLength,
				maxPeptideLength,
				maxMissedCleavages,
				minCharge,
				maxCharge,
				numberOfReportedPeaks,
				addDecoysToBackgound,
				dontRunDecoys,
				percolatorThreshold,
				percolatorProteinThreshold,
				percolatorVersionNumber,
				percolatorTrainingSetSize,
				percolatorTrainingSetThreshold,
				percolatorTrainingIterations,
				alpha,
				beta,
				dataAcquisitionType,
				numberOfThreadsUsed,
				targetWindowCenter,
				precursorWindowSize,
				numberOfQuantitativePeaks,
				minNumOfQuantitativePeaks,
				topNTargetsUsed,
				minIntensity,
				quantifyAcrossSamples,
				true,
				requireVariableMods,
				filterPeaklists,
				doNotUseGlobalFDR,
				precursorIsolationRangeFile,
				percolatorModelFile,
				normalizeByTIC,
				subtractBackground,
				enableAdvancedOptions
		);
	}
}
