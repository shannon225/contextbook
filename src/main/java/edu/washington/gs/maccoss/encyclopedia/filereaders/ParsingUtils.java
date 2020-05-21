package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;

public class ParsingUtils {

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

	public static byte getByte(String parameterName, HashMap<String, String> parameters, byte defaultValue) {
		String value=parameters.get(parameterName);
		if (value==null) {
			return defaultValue;
		}
		try {
			return Byte.parseByte(value);
		} catch (NumberFormatException nfe) {
			throw new EncyclopediaException("Error parsing "+parameterName+" from ["+value+"]", nfe);
		}
	}

	public static double getDouble(String parameterName, HashMap<String, String> parameters, double defaultValue) {
		String value=parameters.get(parameterName);
		if (value==null) {
			return defaultValue;
		}
		try {
			return Double.parseDouble(value);
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
