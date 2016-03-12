package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.HashMap;
import java.util.StringTokenizer;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import gnu.trove.map.hash.TCharFloatHashMap;

public class PecanParameterParser {
	public static HashMap<String,String> getDefaultParameters() {
		HashMap<String, String> map=new HashMap<String, String>();
		map.put("-fixed", "C=57.0214635");
		map.put("-frag", "YONLY");
		map.put("-ptol", "10");
		map.put("-ftol", "10");
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
		map.put("-percolatorLocation", "internal");
		map.put("-deconvoluteOverlappingWindows", "false");
		map.put("-precursorWindowSize", "-1");
		map.put("-numberOfThreadsUsed", Integer.toString(Runtime.getRuntime().availableProcessors()));
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
		final MassTolerance precursorTolerance;
		final MassTolerance fragmentTolerance;
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
		final boolean deconvoluteOverlappingWindows;
		final int numberOfThreadsUsed;
		final float targetWindowCenter;
		final float precursorWindowSize;

		TCharFloatHashMap fixedMods=new TCharFloatHashMap();
		String value=parameters.get("-fixed");
		if (value!=null) {
			try {
				StringTokenizer st=new StringTokenizer(value, ",");
				while (st.hasMoreTokens()) {
					String token=st.nextToken();
					char aa=token.charAt(0);
					float mass=Float.parseFloat(token.substring(2)); // +1 for '=' (could actually be any deliminator)
					fixedMods.put(aa, mass);
				}
			} catch (Exception e) {
				throw new EncyclopediaException("Error parsing fixed modifications from ["+value+"]", e);
			}
		} else {
			fixedMods.put('C', 57.0214635f);
		}
		aaConstants=new AminoAcidConstants(fixedMods);
		
		value=parameters.get("-frag");
		if (value==null) {
			fragType=FragmentationType.YONLY;
		} else if ("CID".equals(value)) {
			fragType=FragmentationType.CID;
		} else if ("ETD".equals(value)) {
			fragType=FragmentationType.ETD;
		} else if ("YONLY".equals(value)) {
			fragType=FragmentationType.YONLY;
		} else {
			throw new EncyclopediaException("Error parsing fragmentation type from ["+value+"]");
		}
		
		value=parameters.get("-ptol");
		if (value==null) {
			precursorTolerance=new MassTolerance(10);
		} else {
			try {
				precursorTolerance=new MassTolerance(Double.parseDouble(value));
			} catch (NumberFormatException nfe) {
				throw new EncyclopediaException("Error parsing precursor tolerance from ["+value+"]", nfe);
			}
		}
		
		value=parameters.get("-ftol");
		if (value==null) {
			fragmentTolerance=new MassTolerance(10);
		} else {
			try {
				fragmentTolerance=new MassTolerance(Double.parseDouble(value));
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
		deconvoluteOverlappingWindows=SearchParameterParser.getBoolean("-deconvoluteOverlappingWindows", parameters, false);
		numberOfThreadsUsed=SearchParameterParser.getInteger("-numberOfThreadsUsed", parameters, Runtime.getRuntime().availableProcessors());
		targetWindowCenter=SearchParameterParser.getFloat("-targetWindowCenter", parameters, -1f);
		precursorWindowSize=SearchParameterParser.getFloat("-precursorWindowSize", parameters, -1f);

		value=parameters.get("-percolatorLocation");
		File percolator=null;
		if (value==null||"internal".equalsIgnoreCase(value)) {
			percolator=null;
		} else {
			percolator=new File(value);
			if (!percolator.exists()||!percolator.canExecute()) {
				Logger.errorLine("Could not execute external Percolator from ["+value+"]. Falling back to internal Percolator");
				percolator=null;
			}
		}
		
		return new PecanSearchParameters(aaConstants, fragType, precursorTolerance, fragmentTolerance, enzyme, minPeptideLength, maxPeptideLength, maxMissedCleavages, minCharge, maxCharge, minEluteTime, numberOfReportedPeaks, addDecoysToBackgound, dontRunDecoys, percolatorThreshold, alpha, beta, percolator, deconvoluteOverlappingWindows, numberOfThreadsUsed, targetWindowCenter, precursorWindowSize);
	}
}
