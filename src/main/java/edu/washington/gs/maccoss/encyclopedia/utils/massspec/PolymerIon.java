package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;

import com.google.common.collect.ComparisonChain;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;

public class PolymerIon extends PrecursorIon {
	int n;
	Polymer p;
	
	public PolymerIon(int n, Polymer p) {
		super(p.getName()+" "+n, p.getMass(n));
		this.n = n;
		this.p = p;
	}
	public static ArrayList<PolymerIon> getAllPolymerProducts(Range r) {
		ArrayList<PolymerIon> ions=new ArrayList<>();
		for (Polymer p : Polymer.polymers) {
			ions.addAll(getPolymerProducts(p, r));
		}
		return ions;
	}
	
	/**
	 * @param r
	 * @return
	 */
	public static ArrayList<PolymerIon> getPolymerProducts(Polymer p, Range r) {
		ArrayList<PolymerIon> ions=new ArrayList<>();
		int n=0;
		while (true) {
			n++;

			PolymerIon ion = new PolymerIon(n, p);
			if (ion.getMass()<r.getStart()) continue;
			if (ion.getMass()>r.getStop()) break;
			
			ions.add(ion);
		}
		return ions;
	}
	
	public int getN() {
		return n;
	}
	public Polymer getPolymer() {
		return p;
	}

	@Override
	public int hashCode() {
		// Note that equal objects will always have identical masses (see below)
		return (int)(getMass()*100.0);
	}
	
	@Override
	public int compareTo(PrecursorIon o) {
		if (o == null|| !(o instanceof PolymerIon)) {
			return 1;
		}
		PolymerIon i=(PolymerIon)o;
		int c=tolerance.compareTo(getMass(), i.getMass());
		if (c!=0) return c;

		// Comparison uses exact mass as well as type and index. Natural ordering will be by mass, with ties settled
		// by type (ordered by declaration order), then index. This will be transitive and consistent with equals, as
		// all comparisons are exact.
		return ComparisonChain.start()
				.compare(p.getName(), i.p.getName())
				.compare(getN(), i.getN())
				.result();
	}
	
	private static final MassTolerance tolerance=new MassTolerance(0.1); // 1 ppm is about the accuracy of floats 
}