package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import gnu.trove.map.hash.TCharFloatHashMap;

public class SearchParameterParser {
	public static HashMap<String,String> getDefaultParameters() {
		HashMap<String, String> map=new HashMap<String, String>();
		map.put("-fixed", "C=57.0214635");
		map.put("-frag", "CID");
		map.put("-ptol", "10");
		map.put("-ftol", "10");
		map.put("-enzyme", "trypsin");
		map.put("-percolatorThreshold", "0.01");
		map.put("-percolatorLocation", "internal");
		map.put("-expectedPeakWidth", "25");
		map.put("-deconvoluteOverlappingWindows", "false");
		map.put("-runPhosphoLocalization", "false");
		return map;
	}
	
	public static SearchParameters getDefaultParametersObject() {
		return parseParameters(getDefaultParameters());
	}
	
	public static SearchParameters parseParameters(File defaultParameters, HashMap<String, String> parameters) {
		HashMap<String, String> map=readFile(defaultParameters);
		map.putAll(parameters);
		return parseParameters(map);
	}
	
	public static SearchParameters parseParameters(HashMap<String, String> parameters) {
		final AminoAcidConstants aaConstants=new AminoAcidConstants(new TCharFloatHashMap());
		final FragmentationType fragType;
		final MassTolerance precursorTolerance;
		final MassTolerance fragmentTolerance;
		final DigestionEnzyme enzyme;
		final float percolatorThreshold;
		final boolean deconvoluteOverlappingWindows;
		final int numberOfThreadsUsed;
		final float targetWindowCenter;
		final float expectedPeakWidth;
		final boolean runPhosphoLocalization;
		
		String value=parameters.get("-frag");
		if (value==null) {
			fragType=FragmentationType.CID;
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

		percolatorThreshold=getFloat("-percolatorThreshold", parameters, 0.01f);
		deconvoluteOverlappingWindows=getBoolean("-deconvoluteOverlappingWindows", parameters, false);
		numberOfThreadsUsed=SearchParameterParser.getInteger("-numberOfThreadsUsed", parameters, Runtime.getRuntime().availableProcessors());
		targetWindowCenter=SearchParameterParser.getFloat("-targetWindowCenter", parameters, -1f);
		expectedPeakWidth=SearchParameterParser.getFloat("-expectedPeakWidth", parameters, 25f);
		runPhosphoLocalization=getBoolean("-runPhosphoLocalization", parameters, false);
		
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
		
		return new SearchParameters(aaConstants, fragType, precursorTolerance, fragmentTolerance, enzyme, percolatorThreshold, percolator, deconvoluteOverlappingWindows, numberOfThreadsUsed, expectedPeakWidth, targetWindowCenter, runPhosphoLocalization);
	}

	public static boolean getBoolean(String parameterName, HashMap<String, String> parameters, boolean defaultValue) {
		String value=parameters.get(parameterName);
		if (value==null) {
			return defaultValue;
		}
		if ("false".equalsIgnoreCase(value)) return false;
		if ("true".equalsIgnoreCase(value)) return true;
		if ("f".equalsIgnoreCase(value)) return false;
		if ("t".equalsIgnoreCase(value)) return true;
		throw new EncyclopediaException("Error parsing "+parameterName+" from ["+value+"]");
	}

	public static int getInteger(String parameterName, HashMap<String, String> parameters, int defaultValue) {
		String value=parameters.get(parameterName);
		if (value==null) {
			return defaultValue;
		}
		try {
			return Integer.parseInt(value);
		} catch (NumberFormatException nfe) {
			throw new EncyclopediaException("Error parsing "+parameterName+" from ["+value+"]", nfe);
		}
	}

	public static float getFloat(String parameterName, HashMap<String, String> parameters, float defaultValue) {
		String value=parameters.get(parameterName);
		if (value==null) {
			return defaultValue;
		}
		try {
			return Float.parseFloat(value);
		} catch (NumberFormatException nfe) {
			throw new EncyclopediaException("Error parsing "+parameterName+" from ["+value+"]", nfe);
		}
	}
	
	public static Pair<String, String> parseEntry(String eachline) {
		String first=eachline.substring(0, eachline.indexOf('=')-1);
		String second=eachline.substring(eachline.indexOf('=')+1);
		Pair<String, String> entry=new Pair<String, String>(first, second);
		return entry;
	}
	
	public static HashMap<String, String> readFile(File f) {
		try {
			BufferedReader in=new BufferedReader(new FileReader(f));

			HashMap<String, String> map=new HashMap<String, String>();
			try {
				String eachline;
				while ((eachline=in.readLine())!=null) {
					if (eachline.trim().length()==0) {
						continue;
					}
					Pair<String, String> entry=parseEntry(eachline);
					
					map.put(entry.x, entry.y);
				}
				return map;

			} catch (IOException ioe) {
				throw new EncyclopediaException("Error parsing parameters from ["+f.getAbsolutePath()+"]");
			} finally {
				if (in!=null) {
					try {
						in.close();
					} catch (IOException ioe) {
						ioe.printStackTrace();
					}
				}
			}
		} catch (IOException ioe) {
			throw new EncyclopediaException("Error parsing parameters from ["+f.getAbsolutePath()+"]");
		}
	}
}
