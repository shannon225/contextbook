package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.list.array.TDoubleArrayList;

/**
 * keeps track of ions/ranges that can't be used for localization/quant
 * NOTE, NOT THREAD SAFE!!!
 * @author searleb
 *
 */
public class FragmentIonBlacklist {
	private final ArrayList<BlacklistedIon> blacklist=new ArrayList<>();
	private final MassTolerance tolerance;

	// mutable:
	private double[] masses=new double[0];
	private Range[] blacklistRanges=new Range[0];
	
	public FragmentIonBlacklist(MassTolerance tolerance) {
		this.tolerance=tolerance;
	}
	
	public int size() {
		return blacklist.size();
	}
	
	public void addIonsToBlacklist(FragmentIonBlacklist toAdd) {
		blacklist.addAll(toAdd.blacklist);
		processNewIons();
	}
	
	public void addIonToBlacklist(double target, Range range) {
		blacklist.add(new BlacklistedIon(target, range));
		processNewIons();
	}

	private void processNewIons() {
		Collections.sort(blacklist);
		
		Pair<double[], Range[]> pairs=toArrays(blacklist);
		masses=pairs.x;
		blacklistRanges=pairs.y;
	}
	
	public boolean isBlacklisted(double target, float rtInSec) {
		int[] indices=tolerance.getIndicies(masses, target);
		for (int i : indices) {
			if (blacklistRanges[i].contains(rtInSec)) {
				return true;
			}
		}
		return false;
	}
	
	public boolean isBlacklisted(double target) {
		int[] indices=tolerance.getIndicies(masses, target);
		if (indices.length>0) return true;
		return false;
	}

	private class BlacklistedIon implements Comparable<BlacklistedIon> {
		public final double mass;
		public final Range blacklist;

		public BlacklistedIon(double mass, Range blacklist) {
			this.mass=mass;
			this.blacklist=blacklist;
		}
		
		@Override
		public int compareTo(BlacklistedIon o) {
			if (o==null) return 1;
			if (mass>o.mass) return 1;
			if (mass<o.mass) return -1;
			return 0;
		}
	}
	
	private static Pair<double[], Range[]> toArrays(Collection<BlacklistedIon> peaks) {
		TDoubleArrayList masses=new TDoubleArrayList();
		ArrayList<Range> blacklistRanges=new ArrayList<Range>();
		for (BlacklistedIon peak : peaks) {
			masses.add(peak.mass);
			blacklistRanges.add(peak.blacklist);
		}
		return new Pair<double[], Range[]>(masses.toArray(), blacklistRanges.toArray(new Range[blacklistRanges.size()]));
	}
}
