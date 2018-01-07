package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import gnu.trove.map.hash.TDoubleIntHashMap;

public abstract class AbstractPecanFragmentationModel extends FragmentationModel {

	public AbstractPecanFragmentationModel(FragmentationModel model) {
		super(model.getMasses(), model.getNeutralLosses(), model.getAas());
	}
	
	public abstract PecanLibraryEntry getPecanSpectrum(byte precursorCharge, double[] sortedBinCounterKeys, TDoubleIntHashMap binCounter, Range fragmentationRange, SearchParameters params, boolean isDecoy);
}
