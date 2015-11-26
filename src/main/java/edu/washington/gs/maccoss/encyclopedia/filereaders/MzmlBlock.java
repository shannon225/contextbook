package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;

public class MzmlBlock {
	public static final MzmlBlock POISON_BLOCK=new MzmlBlock(new ArrayList<PrecursorScan>(), new ArrayList<Stripe>());
	
	private final ArrayList<PrecursorScan> precursors=new ArrayList<PrecursorScan>();
	private final ArrayList<Stripe> stripes=new ArrayList<Stripe>();

	public MzmlBlock(ArrayList<PrecursorScan> precursors, ArrayList<Stripe> stripes) {
		this.precursors.addAll(precursors);
		this.stripes.addAll(stripes);
	}
	
	public ArrayList<PrecursorScan> getPrecursors() {
		return precursors;
	}
	public ArrayList<Stripe> getStripes() {
		return stripes;
	}
}
