package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PecanLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import gnu.trove.map.hash.TDoubleIntHashMap;

public abstract class AbstractPecanFragmentationModel extends FragmentationModel {

	public AbstractPecanFragmentationModel(String modifiedSequence, AminoAcidConstants aaConstants) {
		super(modifiedSequence, aaConstants);
	}
	
	public abstract PecanLibraryEntry getPecanSpectrum(byte precursorCharge, double[] sortedBinCounterKeys, TDoubleIntHashMap binCounter, SearchParameters params, boolean isDecoy);
}
