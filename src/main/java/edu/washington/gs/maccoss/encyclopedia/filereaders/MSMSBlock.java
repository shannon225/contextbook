package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.util.ArrayList;
import java.util.List;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;

public class MSMSBlock {
	public static final MSMSBlock POISON_BLOCK=new MSMSBlock(new ArrayList<PrecursorScan>(), new ArrayList<FragmentScan>());

	private final ArrayList<PrecursorScan> precursors=new ArrayList<PrecursorScan>();
	private final ArrayList<FragmentScan> stripes=new ArrayList<FragmentScan>();

	public MSMSBlock(List<PrecursorScan> precursors, List<FragmentScan> stripes) {
		this.precursors.addAll(precursors);
		this.stripes.addAll(stripes);
	}
	
	public ArrayList<PrecursorScan> getPrecursorScans() {
		return precursors;
	}
	public ArrayList<FragmentScan> getFragmentScans() {
		return stripes;
	}
}
