package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.StringTokenizer;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import gnu.trove.map.hash.TCharFloatHashMap;

public class SearchParameterParser {
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
		return map;
	}
	
	public static SearchParameters parseParameters(File defaultParameters, HashMap<String, String> parameters) {
		HashMap<String, String> map=readFile(defaultParameters);
		map.putAll(parameters);
		return parseParameters(map);
	}
	
	public static SearchParameters parseParameters(HashMap<String, String> parameters) {
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

		minPeptideLength=getInteger("-minLength", parameters, 5);
		maxPeptideLength=getInteger("-maxLength", parameters, 100);
		maxMissedCleavages=getInteger("-maxMissedCleavage", parameters, 1);
		minCharge=(byte)getInteger("-minCharge", parameters, 2);
		maxCharge=(byte)getInteger("-maxCharge", parameters, 3);
		minEluteTime=getInteger("-minEluteTime", parameters, 12);
		numberOfReportedPeaks=getInteger("-numberOfReportedPeaks", parameters, 1);
		addDecoysToBackgound=getBoolean("-addDecoysToBackground", parameters, false);
		dontRunDecoys=getBoolean("-dontRunDecoys", parameters, false);
		return new SearchParameters(aaConstants, fragType, precursorTolerance, fragmentTolerance, enzyme, minPeptideLength, maxPeptideLength, maxMissedCleavages, minCharge, maxCharge, minEluteTime, numberOfReportedPeaks, addDecoysToBackgound, dontRunDecoys);
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
	
	public static Pair<String, String> parseEntry(String eachline) {
		String first=eachline.substring(0, eachline.indexOf('=')-1);
		String second=eachline.substring(eachline.indexOf('=')+1);
		Pair<String, String> entry=new Pair<String, String>(first, second);
		return entry;
	}
}
