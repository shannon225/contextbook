package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ChromatogramLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.set.hash.TIntHashSet;

public class IonAccountingTest {
	public static void main(String[] args) throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		StripeFileInterface dia=StripeFileGenerator.getFile(new File("/Users/searleb/Documents/teaching/encyclopedia/quantitative samples/23aug2017_hela_serum_timecourse_wide_1a.mzML"), parameters, true);
		Map<Range, WindowData> ranges=dia.getRanges();
		
		LibraryInterface diaLib=BlibToLibraryConverter.getFile(new File("/Users/searleb/Documents/teaching/encyclopedia/quantitative samples/23aug2017_hela_serum_timecourse_wide_1a.mzML.elib"));
		
		ArrayList<Range> rangeList=new ArrayList<Range>(ranges.keySet());
		Collections.sort(rangeList);
		
		///////////////////////////////////////////////////
		//rangeList=new ArrayList<>();
		//rangeList.add(new Range(440.25, 440.27));
		///////////////////////////////////////////////////
		
		for (Range target : rangeList) {
			System.out.println("Precursor: "+target);
			final ArrayList<LibraryEntry> diaEntries=diaLib.getEntries(target, false, parameters.getAAConstants());
			
			HashMap<String, ArrayList<LibraryEntry>> lastNSeqs=new HashMap<>();
			for (LibraryEntry entry : diaEntries) {
				if (entry.getPeptideSeq().length()>6) {
					String key=entry.getPeptideSeq().substring(entry.getPeptideSeq().length()-6);
					ArrayList<LibraryEntry> list=lastNSeqs.get(key);
					if (list==null) {
						list=new ArrayList<>();
						lastNSeqs.put(key, list);
					}
					list.add(entry);
				}
			}
			
			for (ArrayList<LibraryEntry> group : lastNSeqs.values()) {
				if (group.size()>1) {
					boolean found=false;
					double mass=-1.0;
					for (LibraryEntry entry : group) {
						if (mass==-1.0) {
							mass=entry.getPrecursorMZ();
						} else if (Math.round(mass)!=Math.round(entry.getPrecursorMZ())) {
							found=true;
							break;
						}
					}
					if (found) {
						for (LibraryEntry entry : group) {
							System.out.println(entry.getPrecursorMZ()+", "+entry.getPeptideModSeq());
						}
						System.out.println();
					}
				}
			}
		}
	}
	public static void main2(String[] args) throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		StripeFileInterface dia=StripeFileGenerator.getFile(new File("/Users/searleb/Documents/teaching/encyclopedia/quantitative samples/23aug2017_hela_serum_timecourse_wide_1a.mzML"), parameters, true);
		Map<Range, WindowData> ranges=dia.getRanges();
		
		LibraryInterface diaLib=BlibToLibraryConverter.getFile(new File("/Users/searleb/Documents/teaching/encyclopedia/quantitative samples/23aug2017_hela_serum_timecourse_wide_1a.mzML.elib"));
		
		ArrayList<Range> rangeList=new ArrayList<Range>(ranges.keySet());
		Collections.sort(rangeList);
		
		///////////////////////////////////////////////////
		//rangeList=new ArrayList<>();
		//rangeList.add(new Range(440.25, 440.27));
		///////////////////////////////////////////////////
		
		for (Range target : rangeList) {
			System.out.println("Precursor: "+target);
			
			final ArrayList<LibraryEntry> diaEntries=diaLib.getEntries(target, false, parameters.getAAConstants());
			Collections.sort(diaEntries);
			//System.out.println("num peptides considered: "+diaEntries.size());
			
			IonAccounter accounter=new IonAccounter(new MassTolerance(10.0));
			int numberOfPeptidesWithAnIon=0;
			for (int peptideIndex = 0; peptideIndex < diaEntries.size(); peptideIndex++) {
				LibraryEntry entry=diaEntries.get(peptideIndex);
				if (entry instanceof ChromatogramLibraryEntry) {
					Range range=((ChromatogramLibraryEntry)entry).getRtRange();
					float[] correlation=entry.getCorrelationArray();
					double[] mzs=entry.getMassArray();
					
					int ionCount=0;
					for (int i=0; i<mzs.length; i++) {
						if (mzs[i]>200&&correlation[i]>TransitionRefiner.identificationCorrelationThreshold) {
							accounter.addIonToList(mzs[i], range, peptideIndex);
							ionCount++;
						}
					}
					if (ionCount>=3) {
						numberOfPeptidesWithAnIon++;
					}
				}
			}
			//System.out.println("num peptides with at least three high quality ions: "+numberOfPeptidesWithAnIon);
			
			accounter.runIonAccounting();
	
			int numberOfPeptidesWithAUniqueIon=0;
			TIntArrayList numberOfUniqueIonsPerPeptide=new TIntArrayList();
			TIntArrayList numberOfIonsPerPeptide=new TIntArrayList();
			for (int peptideIndex = 0; peptideIndex < diaEntries.size(); peptideIndex++) {
				LibraryEntry entry=diaEntries.get(peptideIndex);
				if (entry instanceof ChromatogramLibraryEntry) {
					Range range=((ChromatogramLibraryEntry)entry).getRtRange();
					float[] correlation=entry.getCorrelationArray();
					double[] mzs=entry.getMassArray();
					
					TIntIntHashMap conflictingIDs=new TIntIntHashMap(); // no entry value is 0 by default
					int numUniqueIons=0;
					int numIons=0;
					for (int i=0; i<mzs.length; i++) {
						if (mzs[i]>200&&correlation[i]>TransitionRefiner.identificationCorrelationThreshold) {
							numIons++;
							if (!accounter.isIonOverlapping(mzs[i], range.getMiddle(), peptideIndex)) {
								numUniqueIons++;
							} else {
								int[] overlapping=accounter.getOverlappingIDs(mzs[i], range.getMiddle(), peptideIndex);
								for (int j = 0; j < overlapping.length; j++) {
									int original=conflictingIDs.get(overlapping[j]);
									conflictingIDs.put(overlapping[j], original+1);
								}
							}
						}
					}
					numberOfUniqueIonsPerPeptide.add(numUniqueIons);
					numberOfIonsPerPeptide.add(numIons);
					if (numUniqueIons>=3) {
						numberOfPeptidesWithAUniqueIon++;
						
					} else if (numIons>=6) {
						// as least as many ions as AAs, but no unique ions --> check for conflicts
						TIntHashSet bestIDs=new TIntHashSet();
						int bestCount=0;
	
						for (TIntIntIterator it = conflictingIDs.iterator(); it.hasNext();) {
							it.advance();
							if (it.value() > bestCount) {
								bestCount = it.value();
								bestIDs.clear();
							}
							if (it.value() == bestCount) {
								bestIDs.add(it.key());
							}
						}
						if (bestIDs.size()==1) {
							// only one problematic peptide
							LibraryEntry altEntry=diaEntries.get(bestIDs.toArray()[0]);
							if (!entry.getPeptideSeq().replace('I', 'L').equals(altEntry.getPeptideSeq().replace('I', 'L'))&&Math.round(entry.getPrecursorMZ())!=Math.round(altEntry.getPrecursorMZ())) {
								System.out.println(numUniqueIons+", "+numIons+", "+bestCount+":\t"+entry.getPeptideModSeq()+" ("+General.toString(entry.getAccessions())+")"+" <-versus-> "+altEntry.getPeptideModSeq()+" ("+General.toString(altEntry.getAccessions())+")");
							}
						}
					}
				}
			}
			System.out.println("num peptides with three unique ions: "+numberOfPeptidesWithAUniqueIon+" / "+numberOfPeptidesWithAnIon+" / "+diaEntries.size());
			for (int i=0; i<numberOfUniqueIonsPerPeptide.size(); i++) {
				//System.out.println(numberOfUniqueIonsPerPeptide.get(i)+"\t"+numberOfIonsPerPeptide.get(i));
				LibraryEntry entry=diaEntries.get(i);
				if (numberOfUniqueIonsPerPeptide.get(i)==0&&numberOfIonsPerPeptide.get(i)>6) {
					//System.out.println(numberOfUniqueIonsPerPeptide.get(i)+"\t"+numberOfIonsPerPeptide.get(i)+" --> "+entry.getPeptideModSeq()+" ("+General.toString(entry.getAccessions())+")");
				}
			}
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
		
		public int[] getOverlappingIDs(double target, float rtInSec, int peptideID) {
			TIntHashSet otherIDs=new TIntHashSet();
			int[] indices=tolerance.getIndicies(masses, target);
			for (int i : indices) {
				if (blacklistRanges[i].contains(rtInSec)) {
					if (peptideIDs[i]!=peptideID) {
						otherIDs.add(peptideIDs[i]);
					}
				}
			}
			return otherIDs.toArray();
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
