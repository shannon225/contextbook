package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaToPrositCSVParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;

import java.util.HashMap;

public class FastaToPrositCSVParametersParser {

	public static HashMap<String,String> getDefaultParameters() {
		FastaToPrositCSVParameters defaultParameters = parseParameters(new HashMap<>());
		return defaultParameters.toParameterMap();
	}

	public static FastaToPrositCSVParameters parseParameters(HashMap<String, String> parameters) {
		final int defaultNCE;
		final byte defaultCharge;
		final int minCharge;
		final int maxCharge;
		final int maxMissedCleavage;
		final double minMz;
		final double maxMz;
		final DigestionEnzyme enzyme;

		defaultNCE= FastaToPrositCSVParametersParser.getInteger("-defaultNCE", parameters, FastaToPrositCSVParameters.DEFAULT_DEFAULT_NCE);
		defaultCharge= FastaToPrositCSVParametersParser.getByte("-defaultCharge", parameters, FastaToPrositCSVParameters.DEFAULT_DEFAULT_CHARGE);
		minCharge= FastaToPrositCSVParametersParser.getInteger("-minCharge", parameters, FastaToPrositCSVParameters.DEFAULT_MIN_CHARGE);
		maxCharge= FastaToPrositCSVParametersParser.getInteger("-maxCharge", parameters, FastaToPrositCSVParameters.DEFAULT_MAX_CHARGE);
		maxMissedCleavage= FastaToPrositCSVParametersParser.getInteger("-maxMissedCleavage", parameters, FastaToPrositCSVParameters.DEFAULT_MAX_MISSED_CLEAVAGE);
		minMz= FastaToPrositCSVParametersParser.getDouble("-minMz", parameters, FastaToPrositCSVParameters.DEFAULT_MIN_MZ);
		maxMz= FastaToPrositCSVParametersParser.getDouble("-maxMz", parameters, FastaToPrositCSVParameters.DEFAULT_MAX_MZ);

		String value=parameters.get("-enzyme");
		if (value==null) {
			enzyme= FastaToPrositCSVParameters.DEFAULT_ENZYME;
		} else {
			enzyme= DigestionEnzyme.getEnzyme(value);
		}

		return new FastaToPrositCSVParameters(defaultNCE, defaultCharge, minCharge, maxCharge, maxMissedCleavage, minMz, maxMz, enzyme);
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
		} else {
			try {
				return Double.parseDouble(value);
			} catch (NumberFormatException nfe) {
				throw new EncyclopediaException("Error parsing fragment tolerance from ["+value+"]", nfe);
			}
		}
	}
}
