package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.NoSuchElementException;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Iterables;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ScoringBreadthType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
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
		map.put("-percolatorThreshold", "0.01");
		map.put("-percolatorProteinThreshold", "0.01");
		map.put("-percolatorVersionNumber", Byte.toString(PercolatorExecutor.DEFAULT_VERSION_NUMBER));
		map.put("-expectedPeakWidth", "25");
		map.put("-acquisition", DataAcquisitionType.toName(DataAcquisitionType.DIA));
		map.put("-localizationModification", PeptideModification.NO_MODIFICATION_NAME);
		map.put("-scoringBreadthType", ScoringBreadthType.ENTIRE_RT_WINDOW.toShortname());
		map.put("-numberOfExtraDecoyLibrariesSearched", "0.0");
		map.put("-numberOfQuantitativePeaks", "5");
		map.put("-minNumOfQuantitativePeaks", "3");
		map.put("-minQuantitativeIonNumber", "3");
		map.put("-verifyModificationIons", "true");
		map.put("-minIntensity", "-1.0");
		map.put("-rtWindowInMin", "-1.0");
        map.put("-filterPeaklists", "false");
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
		map.put("-numberOfQuantitativePeaks", "5");
		map.put("-minNumOfQuantitativePeaks", "3");
		map.put("-minQuantitativeIonNumber", "3");
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
		final DataAcquisitionType dataAcquisitionType;
		final int numberOfThreadsUsed;
		final float targetWindowCenter;
		final float expectedPeakWidth;
		final float precursorWindowSize;
		final int numberOfQuantitativePeaks;
		final int minNumOfQuantitativePeaks;
		final float minIntensity;
		final float numberOfExtraDecoyLibrariesSearched;
		final int percolatorVersionNumber;
		final Optional<PeptideModification> localizationModification;
		final ScoringBreadthType breadthType;
		final boolean quantifyAcrossSamples;
		final boolean verifyModificationIons;
		final float rtWindowInMin;
        final boolean filterPeaklists;
        final boolean doNotUseGlobalFDR;
		
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

		percolatorThreshold=getFloat("-percolatorThreshold", parameters, 0.01f);
		percolatorProteinThreshold=getFloat("-percolatorProteinThreshold", parameters, 0.01f);
		numberOfThreadsUsed=SearchParameterParser.getInteger("-numberOfThreadsUsed", parameters, Runtime.getRuntime().availableProcessors());
		targetWindowCenter=SearchParameterParser.getFloat("-targetWindowCenter", parameters, -1f);
		precursorWindowSize=SearchParameterParser.getFloat("-precursorWindowSize", parameters, -1f);
		expectedPeakWidth=SearchParameterParser.getFloat("-expectedPeakWidth", parameters, 25f);
		numberOfQuantitativePeaks=SearchParameterParser.getInteger("-numberOfQuantitativePeaks", parameters, 5);
		minNumOfQuantitativePeaks=SearchParameterParser.getInteger("-minNumOfQuantitativePeaks", parameters, 3);
		minIntensity=SearchParameterParser.getFloat("-minIntensity", parameters, -1.0f);
		rtWindowInMin=SearchParameterParser.getFloat("-rtWindowInMin", parameters, -1f);
		
		percolatorVersionNumber=SearchParameterParser.getInteger("-percolatorVersionNumber", parameters, 3);
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
				Logger.errorLine("Falling back to scoring breadth type: "+ScoringBreadthType.ENTIRE_RT_WINDOW.toShortname());
				type=ScoringBreadthType.ENTIRE_RT_WINDOW;
			}
			breadthType=type;
		} else {
			breadthType=ScoringBreadthType.ENTIRE_RT_WINDOW;
		}
		
		float tempNumberOfExtraDecoyLibrariesSearched=SearchParameterParser.getFloat("-numberOfExtraDecoyLibrariesSearched", parameters, 0.0f);
		if (tempNumberOfExtraDecoyLibrariesSearched<0.0f) {
			Logger.errorLine("-numberOfExtraDecoyLibrariesSearched cannot be less than 0%! Using 0% extra decoys.");
			numberOfExtraDecoyLibrariesSearched=0.0f;
		} else {
			numberOfExtraDecoyLibrariesSearched=tempNumberOfExtraDecoyLibrariesSearched;
		}
		quantifyAcrossSamples=SearchParameterParser.getBoolean("-quantifyAcrossSamples", parameters, false);
		verifyModificationIons=SearchParameterParser.getBoolean("-verifyModificationIons", parameters, true);
        filterPeaklists=SearchParameterParser.getBoolean("-filterPeaklists", parameters, false);
        doNotUseGlobalFDR=SearchParameterParser.getBoolean("-doNotUseGlobalFDR", parameters, false);

		return new SearchParameters(aaConstants, fragType, precursorTolerance, precursorOffsetPPM, precursorIsolationMargin, fragmentTolerance, fragmentOffsetPPM, libraryFragmentTolerance, enzyme, percolatorThreshold, percolatorProteinThreshold, percolatorVersionNumber, dataAcquisitionType, numberOfThreadsUsed, expectedPeakWidth,
				targetWindowCenter, precursorWindowSize, numberOfQuantitativePeaks, minNumOfQuantitativePeaks, minIntensity, localizationModification, breadthType, numberOfExtraDecoyLibrariesSearched, quantifyAcrossSamples, verifyModificationIons, rtWindowInMin, filterPeaklists, doNotUseGlobalFDR);
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
