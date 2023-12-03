package edu.washington.gs.maccoss.encyclopedia.datastructures.parameters;

import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassErrorUnitType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

public enum InstrumentSpecificSearchParameters {
	IontrapIontrap, OrbitrapIontrap, OrbitrapOrbitrap, ToF, OrbitrapAstral;
	
	public static final InstrumentSpecificSearchParameters[] INSTRUMENTS={OrbitrapOrbitrap, ToF, OrbitrapAstral, IontrapIontrap, OrbitrapIontrap};
	
	public static HashMap<String, String> checkParameters(HashMap<String, String> params) {
		String val=params.get(SearchParameters.INSTRUMENT);
		if (val==null) return params;
		
		InstrumentSpecificSearchParameters instrument=fromString(val);
		return instrument.underwriteParameters(params);
	}

	/**
	 * prefer parameters previously set
	 * @param params
	 * @return
	 */
	public HashMap<String, String> underwriteParameters(HashMap<String, String> params) {
		HashMap<String, String> instrumentSpecific=overwriteParameters(new HashMap<>());
		instrumentSpecific.putAll(params);
		return instrumentSpecific;
	}

	public HashMap<String, String> overwriteParameters(SearchParameters parameters) {
		HashMap<String, String> params=parameters.toParameterMap();
		return overwriteParameters(params);
	}
	
	/**
	 * prefer instrument-specific parameters
	 * @param params
	 * @return
	 */
	public HashMap<String, String> overwriteParameters(HashMap<String, String> params) {
		
		switch (this) {
			case OrbitrapOrbitrap:
				params.put("-ptol", "10");
				params.put("-ftol", "10");
				params.put("-ptolunits", MassErrorUnitType.PPM.toString());
				params.put("-ftolunits", MassErrorUnitType.PPM.toString());
				params.put(SearchParameters.SUBTRACT_BACKGROUND, Boolean.TRUE.toString());
				params.put("-minIntensity", "-1");
				params.put("-filterPeaklists", Boolean.FALSE.toString());
				break;
				
			case ToF:
				params.put("-ptol", "10");
				params.put("-ftol", "10");
				params.put("-ptolunits", MassErrorUnitType.PPM.toString());
				params.put("-ftolunits", MassErrorUnitType.PPM.toString());
				params.put(SearchParameters.SUBTRACT_BACKGROUND, Boolean.TRUE.toString());
				params.put("-minIntensity", "125");
				params.put("-filterPeaklists", Boolean.TRUE.toString());
				break;
				
			case OrbitrapAstral:
				params.put("-ptol", "10");
				params.put("-ftol", "10");
				params.put("-ptolunits", MassErrorUnitType.PPM.toString());
				params.put("-ftolunits", MassErrorUnitType.PPM.toString());
				params.put(SearchParameters.SUBTRACT_BACKGROUND, Boolean.TRUE.toString());
				params.put("-minIntensity", "500");
				params.put("-filterPeaklists", Boolean.TRUE.toString());
				break;
				
			case IontrapIontrap:
				params.put("-ptol", "0.4");
				params.put("-ftol", "0.4");
				params.put("-ptolunits", MassErrorUnitType.AMU.toString());
				params.put("-ftolunits", MassErrorUnitType.AMU.toString());
				params.put(SearchParameters.SUBTRACT_BACKGROUND, Boolean.TRUE.toString());
				params.put("-minIntensity", "-1");
				params.put("-filterPeaklists", Boolean.TRUE.toString());
				break;

			case OrbitrapIontrap:
				params.put("-ptol", "10");
				params.put("-ftol", "0.4");
				params.put("-ptolunits", MassErrorUnitType.PPM.toString());
				params.put("-ftolunits", MassErrorUnitType.AMU.toString());
				params.put(SearchParameters.SUBTRACT_BACKGROUND, Boolean.TRUE.toString());
				params.put("-minIntensity", "-1");
				params.put("-filterPeaklists", Boolean.TRUE.toString());
				break;
	
			default:
				Logger.errorLine("Unknown parameters for instrument type: "+this);
				break;
		}
		
		return params;
	}
	
	public MassTolerance getPrecursorTolerance() {
		HashMap<String, String> params=overwriteParameters(new HashMap<>());
		MassErrorUnitType toleranceType=MassErrorUnitType.getUnitType(params.get("-ptolunits"));
		return new MassTolerance(Double.parseDouble(params.get("-ptol")), toleranceType);
	}
	
	public MassTolerance getFragmentTolerance() {
		HashMap<String, String> params=overwriteParameters(new HashMap<>());
		MassErrorUnitType toleranceType=MassErrorUnitType.getUnitType(params.get("-ftolunits"));
		return new MassTolerance(Double.parseDouble(params.get("-ftol")), toleranceType);
	}
	
	public String toString() {
		switch (this) {
		case OrbitrapOrbitrap:
			return "Orbitrap/Orbitrap";
		case ToF:
			return "ToF";
		case OrbitrapAstral:
			return "Astral";
		case IontrapIontrap:
			return "Iontrap/Iontrap";
		case OrbitrapIontrap:
			return "Orbitrap/Iontrap";
		default:
			return "Unknown";
		}
	}
	
	public static InstrumentSpecificSearchParameters fromString(String s) {
		s=s.toLowerCase();
		
		if (s.equals("OrbiTrap/OrbiTrap".toLowerCase())) {
			return OrbitrapOrbitrap;
		} else if (s.equals("OrbiTrap".toLowerCase())) {
			return OrbitrapOrbitrap;
		} else if (s.equals("Orbi".toLowerCase())) {
			return OrbitrapOrbitrap;
		} else if (s.equals("OT".toLowerCase())) {
			return OrbitrapOrbitrap;
		} else if (s.equals("OrbiTrap/Astral".toLowerCase())) {
			return OrbitrapAstral;
		} else if (s.equals("Astral".toLowerCase())) {
			return OrbitrapAstral;
		} else if (s.equals("IonTrap/IonTrap".toLowerCase())) {
			return IontrapIontrap;
		} else if (s.equals("IonTrap".toLowerCase())) {
			return IontrapIontrap;
		} else if (s.equals("LIT".toLowerCase())) {
			return IontrapIontrap;
		} else if (s.equals("IT".toLowerCase())) {
			return IontrapIontrap;
		} else if (s.equals("ToF".toLowerCase())) {
			return ToF;
		} else if (s.equals("OrbitrapIontrap".toLowerCase())) {
			return OrbitrapIontrap;
		} else if (s.equals("OT/IT".toLowerCase())) {
			return OrbitrapIontrap;
		}

		throw new EncyclopediaException("Unknown instrument type: "+s);
	}
}
