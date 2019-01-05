package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;

public class MSMSBlock {
	public static final MSMSBlock POISON_BLOCK=new MSMSBlock(new ArrayList<PrecursorScan>(), new ArrayList<FragmentScan>());
	
	private final ArrayList<PrecursorScan> precursors=new ArrayList<PrecursorScan>();
	private final ArrayList<FragmentScan> stripes=new ArrayList<FragmentScan>();

	public MSMSBlock(ArrayList<PrecursorScan> precursors, ArrayList<FragmentScan> stripes) {
		this.precursors.addAll(precursors);
		this.stripes.addAll(stripes);
	}
	
	public ArrayList<PrecursorScan> getPrecursors() {
		return precursors;
	}
	public ArrayList<FragmentScan> getStripes() {
		return stripes;
	}
}
