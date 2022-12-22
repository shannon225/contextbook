package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaToPrositCSVParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;

public class FastaToPrositCSVParametersParser {

	public static HashMap<String,String> getDefaultParameters() {
		FastaToPrositCSVParameters defaultParameters = parseParameters(new HashMap<>());
		return defaultParameters.toParameterMap();
	}

	public static FastaToPrositCSVParameters parseParameters(HashMap<String, String> parameters) {
		final int defaultNCE;
		final byte defaultCharge;
		final byte minCharge;
		final byte maxCharge;
		final int maxMissedCleavage;
		final double minMz;
		final double maxMz;
		final DigestionEnzyme enzyme;
		final boolean adjustNCEForDIA;
		final boolean addDecoys;

		defaultNCE= ParsingUtils.getInteger("-defaultNCE", parameters, FastaToPrositCSVParameters.DEFAULT_DEFAULT_NCE);
		defaultCharge= ParsingUtils.getByte("-defaultCharge", parameters, FastaToPrositCSVParameters.DEFAULT_DEFAULT_CHARGE);
		minCharge= ParsingUtils.getByte("-minCharge", parameters, FastaToPrositCSVParameters.DEFAULT_MIN_CHARGE);
		maxCharge= ParsingUtils.getByte("-maxCharge", parameters, FastaToPrositCSVParameters.DEFAULT_MAX_CHARGE);
		maxMissedCleavage= ParsingUtils.getInteger("-maxMissedCleavage", parameters, FastaToPrositCSVParameters.DEFAULT_MAX_MISSED_CLEAVAGE);
		minMz= ParsingUtils.getDouble("-minMz", parameters, FastaToPrositCSVParameters.DEFAULT_MIN_MZ);
		maxMz= ParsingUtils.getDouble("-maxMz", parameters, FastaToPrositCSVParameters.DEFAULT_MAX_MZ);
		adjustNCEForDIA= ParsingUtils.getBoolean("-adjustNCEForDIA", parameters, FastaToPrositCSVParameters.DEFAULT_ADJUST_NCE_FOR_DIA);
		addDecoys= ParsingUtils.getBoolean("-addDecoys", parameters, FastaToPrositCSVParameters.DEFAULT_ADD_DECOYS);

		String value=parameters.get("-enzyme");
		if (value==null) {
			enzyme= FastaToPrositCSVParameters.DEFAULT_ENZYME;
		} else {
			enzyme= DigestionEnzyme.getEnzyme(value);
		}

		return new FastaToPrositCSVParameters(defaultNCE, defaultCharge, minCharge, maxCharge, maxMissedCleavage, minMz, maxMz, enzyme, adjustNCEForDIA, addDecoys);
	}
}
