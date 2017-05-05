package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.HashMap;
import java.util.StringTokenizer;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
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
		map.put("-enzyme", "trypsin");
		map.put("-minLength", "5");
		map.put("-maxLength", "100");
		map.put("-maxMissedCleavage", "1");
		map.put("-minCharge", "2");
		map.put("-maxCharge", "3");
		map.put("-minEluteTime", "12");
		map.put("-numberOfReportedPeaks", "1");
		map.put("-addDecoysToBackground", "false");
		map.put("-dontRunDecoys", "false");
		map.put("-percolatorThreshold", "0.01");
		map.put("-alpha", "1.8");
		map.put("-beta", "0.4");
		map.put("-percolatorVersionNumber", Byte.toString(PercolatorExecutor.DEFAULT_VERSION_NUMBER));
		map.put("-acquisition", "overlapping dia");
		map.put("-precursorWindowSize", "-1");
		map.put("-numberOfThreadsUsed", Integer.toString(Runtime.getRuntime().availableProcessors()));
		map.put("-numberOfQuantitativePeaks", "5");
		return map;
	}
	
	public static SearchParameters getDefaultParametersObject() {
		return parseParameters(getDefaultParameters());
	}
	
	public static SearchParameters parseParameters(File defaultParameters, HashMap<String, String> parameters) {
		HashMap<String, String> map=SearchParameterParser.readFile(defaultParameters);
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
		final DigestionEnzyme enzyme;
		final int minPeptideLength;
		final int maxPeptideLength;
		final int maxMissedCleavages;
		final byte minCharge;
		final byte maxCharge;
		final int minEluteTime;
		final int numberOfReportedPeaks;
		final boolean addDecoysToBackgound;
		final boolean dontRunDecoys;
		final float percolatorThreshold;
		final float alpha;
		final float beta;
		final DataAcquisitionType dataAcquisitionType;
		final int numberOfThreadsUsed;
		final float targetWindowCenter;
		final float precursorWindowSize;
		final int numberOfQuantitativePeaks;
		final int percolatorVersionNumber;

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
			fragType=FragmentationType.YONLY;
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
		
		value=parameters.get("-enzyme");
		if (value==null) {
			enzyme=DigestionEnzyme.getEnzyme("trypsin");
		} else {
			enzyme=DigestionEnzyme.getEnzyme(value);
		}

		minPeptideLength=SearchParameterParser.getInteger("-minLength", parameters, 5);
		maxPeptideLength=SearchParameterParser.getInteger("-maxLength", parameters, 100);
		maxMissedCleavages=SearchParameterParser.getInteger("-maxMissedCleavage", parameters, 1);
		minCharge=(byte)SearchParameterParser.getInteger("-minCharge", parameters, 2);
		maxCharge=(byte)SearchParameterParser.getInteger("-maxCharge", parameters, 3);
		minEluteTime=SearchParameterParser.getInteger("-minEluteTime", parameters, 12);
		numberOfReportedPeaks=SearchParameterParser.getInteger("-numberOfReportedPeaks", parameters, 1);
		addDecoysToBackgound=SearchParameterParser.getBoolean("-addDecoysToBackground", parameters, false);
		dontRunDecoys=SearchParameterParser.getBoolean("-dontRunDecoys", parameters, false);
		percolatorThreshold=SearchParameterParser.getFloat("-percolatorThreshold", parameters, 0.01f);
		alpha=SearchParameterParser.getFloat("-alpha", parameters, 1.8f);
		beta=SearchParameterParser.getFloat("-beta", parameters, 0.4f);
		numberOfThreadsUsed=SearchParameterParser.getInteger("-numberOfThreadsUsed", parameters, Runtime.getRuntime().availableProcessors());
		targetWindowCenter=SearchParameterParser.getFloat("-targetWindowCenter", parameters, -1f);
		precursorWindowSize=SearchParameterParser.getFloat("-precursorWindowSize", parameters, -1f);
		numberOfQuantitativePeaks=SearchParameterParser.getInteger("-numberOfQuantitativePeaks", parameters, 5);
		percolatorVersionNumber=SearchParameterParser.getInteger("-percolatorVersionNumber", parameters, 3);
		
		return new PecanSearchParameters(aaConstants, fragType, precursorTolerance, precursorOffsetPPM, fragmentTolerance, fragmentOffsetPPM, enzyme, minPeptideLength, maxPeptideLength, maxMissedCleavages, minCharge, maxCharge, minEluteTime, numberOfReportedPeaks, addDecoysToBackgound, dontRunDecoys, percolatorThreshold, alpha, beta, percolatorVersionNumber, dataAcquisitionType, numberOfThreadsUsed, targetWindowCenter, precursorWindowSize, numberOfQuantitativePeaks);
	}
}
