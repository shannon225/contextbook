package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ChromatogramLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;

public class IonAccountingTest {
	public static void main(String[] args) throws Exception {
		LibraryInterface diaLib=BlibToLibraryConverter.getFile(new File("/Users/searleb/Documents/teaching/encyclopedia/quantitative samples/23aug2017_hela_serum_timecourse_wide_1a.mzML.elib"));
		final ArrayList<LibraryEntry> diaEntries=diaLib.getAllEntries(false, new AminoAcidConstants());
		System.out.println("num peptides considered: "+diaEntries.size());
		
		IonAccounter accounter=new IonAccounter(new MassTolerance(10.0));
		int peptideID=0;
		int numberOfPeptidesWithAnIon=0;
		for (LibraryEntry entry : diaEntries) {
			if (entry instanceof ChromatogramLibraryEntry) {
				peptideID=peptideID+1;
				
				Range range=((ChromatogramLibraryEntry)entry).getRtRange();
				float[] correlation=entry.getCorrelationArray();
				double[] mzs=entry.getMassArray();
				
				boolean added=false;
				for (int i=0; i<mzs.length; i++) {
					if (mzs[i]>200&&correlation[i]>TransitionRefiner.identificationCorrelationThreshold) {
						accounter.addIonToList(mzs[i], range, peptideID);
						added=true;
					}
				}
				if (added) {
					numberOfPeptidesWithAnIon++;
				}
			}
		}
		System.out.println("num peptides with at least one high quality ion: "+numberOfPeptidesWithAnIon);
		
		accounter.runIonAccounting();

		peptideID=0;
		int numberOfPeptidesWithAUniqueIon=0;
		TIntArrayList numberOfUniqueIonsPerPeptide=new TIntArrayList();
		TIntArrayList numberOfIonsPerPeptide=new TIntArrayList();
		for (LibraryEntry entry : diaEntries) {
			if (entry instanceof ChromatogramLibraryEntry) {
				peptideID=peptideID+1;
				
				Range range=((ChromatogramLibraryEntry)entry).getRtRange();
				float[] correlation=entry.getCorrelationArray();
				double[] mzs=entry.getMassArray();
				
				int numUniqueIons=0;
				int numIons=0;
				for (int i=0; i<mzs.length; i++) {
					if (mzs[i]>200&&correlation[i]>TransitionRefiner.identificationCorrelationThreshold) {
						numIons++;
						if (!accounter.isIonOverlapping(mzs[i], range.getMiddle(), peptideID)) {
							numUniqueIons++;
						}
					}
				}
				numberOfUniqueIonsPerPeptide.add(numUniqueIons);
				numberOfIonsPerPeptide.add(numIons);
				if (numUniqueIons>0) {
					numberOfPeptidesWithAUniqueIon++;
				}
			}
		}
		System.out.println("num peptides with a unique ion: "+numberOfPeptidesWithAUniqueIon);
		for (int i=0; i<numberOfUniqueIonsPerPeptide.size(); i++) {
			System.out.println(numberOfUniqueIonsPerPeptide.get(i)+"\t"+numberOfIonsPerPeptide.get(i));
			
		}
	}
	
	public static class IonAccounter {
		private final ArrayList<AnnotatedIon> blacklist=new ArrayList<>();
		private final MassTolerance tolerance;

		// mutable:
		private double[] masses=new double[0];
		private Range[] blacklistRanges=new Range[0];
		private int[] peptideIDs=new int[0];
		
		public IonAccounter(MassTolerance tolerance) {
			this.tolerance=tolerance;
		}
		
		public int size() {
			return blacklist.size();
		}
		
		public void addIonToList(IonAccounter toAdd) {
			blacklist.addAll(toAdd.blacklist);
		}
		
		public void addIonToList(double target, Range range, int peptideID) {
			blacklist.add(new AnnotatedIon(target, range, peptideID));
		}

		public void runIonAccounting() {
			Collections.sort(blacklist);
			
			Triplet<double[], Range[], int[]> pairs=toArrays(blacklist);
			masses=pairs.x;
			blacklistRanges=pairs.y;
			peptideIDs=pairs.z;
		}
		
		public boolean isIonOverlapping(double target, float rtInSec, int peptideID) {
			int[] indices=tolerance.getIndicies(masses, target);
			for (int i : indices) {
				if (blacklistRanges[i].contains(rtInSec)) {
					if (peptideIDs[i]!=peptideID) {
						return true;
					}
				}
			}
			return false;
		}
		
		public boolean isBlacklisted(double target) {
			int[] indices=tolerance.getIndicies(masses, target);
			if (indices.length>0) return true;
			return false;
		}

		class AnnotatedIon implements Comparable<AnnotatedIon> {
			public final double mass;
			public final Range blacklist;
			public final int id;

			public AnnotatedIon(double mass, Range blacklist, int id) {
				this.mass=mass;
				this.blacklist=blacklist;
				this.id=id;
			}
			
			@Override
			public int compareTo(AnnotatedIon o) {
				if (o==null) return 1;
				if (mass>o.mass) return 1;
				if (mass<o.mass) return -1;
				if (id>o.id) return 1;
				if (id<o.id) return -1;
				return 0;
			}
		}
		
		private Triplet<double[], Range[], int[]> toArrays(Collection<AnnotatedIon> peaks) {
			TDoubleArrayList masses=new TDoubleArrayList();
			ArrayList<Range> blacklistRanges=new ArrayList<Range>();
			TIntArrayList ids=new TIntArrayList();
			for (AnnotatedIon peak : peaks) {
				masses.add(peak.mass);
				blacklistRanges.add(peak.blacklist);
				ids.add(peak.id);
			}
			return new Triplet<double[], Range[], int[]>(masses.toArray(), blacklistRanges.toArray(new Range[blacklistRanges.size()]), ids.toArray());
		}
	}
}
