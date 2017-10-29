package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.AlleleVariant;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.ExtendedFastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.set.hash.TCharHashSet;

public class DigestionEnzyme {
	private static final char[] AAs="ACDEFGHIKLMNPQRSTVWY".toCharArray();
	private final String name;
	private final String percolatorName;
	private final TCharHashSet nterm;
	private final TCharHashSet cterm;
	
	public static HashSet<DigestionEnzyme> getAvailableEnzymes() {
		HashSet<DigestionEnzyme> enzymes=new HashSet<DigestionEnzyme>();
		enzymes.add(getEnzyme("Trypsin"));
		enzymes.add(getEnzyme("Trypsin/p"));
		enzymes.add(getEnzyme("Lys-C"));
		enzymes.add(getEnzyme("Lys-N"));
		enzymes.add(getEnzyme("Arg-C"));
		enzymes.add(getEnzyme("Chymotrypsin"));
		enzymes.add(getEnzyme("Pepsin A"));
		enzymes.add(getEnzyme("Elastase"));
		enzymes.add(getEnzyme("Thermolysin"));
		enzymes.add(getEnzyme("No Enzyme"));
		return enzymes;
	}
	
	public static DigestionEnzyme getEnzyme(String enzymeName) {
		TCharHashSet n=new TCharHashSet();
		TCharHashSet c=new TCharHashSet();
		if ("Trypsin".equalsIgnoreCase(enzymeName)) {
			n.add('K');
			n.add('R');
			c.addAll(AAs);
			c.remove('P');
			
			return new DigestionEnzyme("Trypsin", "trypsin", n, c);
			
		} else if ("Trypsin/p".equalsIgnoreCase(enzymeName)) {
			n.add('K');
			n.add('R');
			c.addAll(AAs);
			
			return new DigestionEnzyme("Trypsin/p", "trypsinp", n, c);
			
		} else if ("No Enzyme".equalsIgnoreCase(enzymeName)) {
			
			return new DigestionEnzyme("No Enzyme", "no_enzyme", n, c);
			
		} else if ("None".equalsIgnoreCase(enzymeName)) {
			
			return new DigestionEnzyme("No Enzyme", "no_enzyme", n, c);
			
		} else if ("Lys-C".equalsIgnoreCase(enzymeName)) {
			n.add('K');
			c.addAll(AAs);
			c.remove('P');
			
			return new DigestionEnzyme("Lys-C", "lys-c", n, c);
			
		} else if ("Lys-N".equalsIgnoreCase(enzymeName)) {
			n.addAll(AAs);
			c.add('K');
			
			return new DigestionEnzyme("Lys-N", "lys-n", n, c);
			
		} else if ("Arg-C".equalsIgnoreCase(enzymeName)) {
			n.add('R');
			c.addAll(AAs);
			c.remove('P');
			
			return new DigestionEnzyme("Arg-C", "arg-c", n, c);
			
		} else if ("Chymotrypsin".equalsIgnoreCase(enzymeName)) {
			n.add('F');
			n.add('Y');
			n.add('W');
			n.add('L');
			c.addAll(AAs);
			c.remove('P');
			
			return new DigestionEnzyme("Chymotrypsin", "chymotrypsin", n, c);
			
		} else if ("Elastase".equalsIgnoreCase(enzymeName)) {
			n.add('A');
			n.add('V');
			c.addAll(AAs);
			
			return new DigestionEnzyme("Elastase", "elastase", n, c);
			
		} else if ("Thermolysin".equalsIgnoreCase(enzymeName)) {
			c.add('A');
			c.add('F');
			c.add('I');
			c.add('L');
			c.add('M');
			c.add('V');
			n.addAll(AAs);
			n.remove('D');
			n.remove('E');
			
			return new DigestionEnzyme("Thermolysin", "thermolysin", n, c);
			
		} else if ("Pepsin A".equalsIgnoreCase(enzymeName)) {
			n.add('F');
			n.add('L');
			c.addAll(AAs);
			
			return new DigestionEnzyme("Pepsin A", "pepsin", n, c);
		}
		
		throw new EncyclopediaException("Unknown digestion enzyme ["+enzymeName+"]");
	}
	
	DigestionEnzyme(String name, String percolatorName, TCharHashSet nterm, TCharHashSet cterm) {
		this.name=name;
		this.percolatorName=percolatorName;
		this.nterm=nterm;
		this.cterm=cterm;
	}
	
	public String getName() {
		return name;
	}
	
	public String getPercolatorName() {
		return percolatorName;
	}
	
	public boolean isCutSite(char pre, char post) {
		if (nterm.contains(pre)&&cterm.contains(post)) {
			return true;
		}
		return false;
	}
	
	/**
	 * assumes key sequence signatures allow fewer than 10 AAs
	 * @param pre
	 * @return
	 */
	public boolean isTargetPreSite(char pre) {
		if (nterm.size()>=10) return false; 
		if (nterm.contains(pre)) {
			return true;
		}
		return false;
	}

	/**
	 * assumes key sequence signatures allow fewer than 10 AAs
	 * @param post
	 * @return
	 */
	public boolean isTargetPostSite(char post) {
		if (cterm.size()>=10) return false; 
		if (cterm.contains(post)) {
			return true;
		}
		return false;
	}
	
	public String reverseProtein(String sequence) {
		StringBuilder sb=new StringBuilder();
		
		int start=0;
		int stop;

		while (start<sequence.length()) {
			stop=start;
			while ((stop<sequence.length()-1)&&!isCutSite(sequence.charAt(stop), sequence.charAt(stop+1))) {
				stop++;
			}
			String peptide=sequence.substring(start, stop+1);
			sb.append(PeptideUtils.reverse(peptide, this));
			start=stop+1;
		}
		return sb.toString();
	}
	
	//@MoMo 
	public ArrayList<String> digestProtein(FastaEntryInterface entry, int minLength, int maxLength, int maxMissedCleavages, AminoAcidConstants constants) {
		if (entry instanceof ExtendedFastaEntry) {
			return digestProtein(entry.getSequence(), minLength, maxLength, maxMissedCleavages, constants,((ExtendedFastaEntry)entry).getPotentialVariant());
		} else {
			return (digestProtein(entry.getSequence(), minLength, maxLength, maxMissedCleavages, constants,new ArrayList<AlleleVariant>()));
		}
	}

	//@MoMo modified 
	private ArrayList<String> digestProtein(String sequence, int minLength, int maxLength, int maxMissedCleavages, AminoAcidConstants constants, ArrayList<AlleleVariant> variants) {
		int totalAllowedStarts=maxMissedCleavages+1;
		ArrayList<String> peptides=new ArrayList<String>();
		TIntArrayList starts=new TIntArrayList();
		starts.add(0);
		int stop;

		// digestion for canonical sequence
		while (starts.get(starts.size()-1)<sequence.length()) {
			stop=starts.get(starts.size()-1);
			while ((stop<sequence.length()-1)&&!isCutSite(sequence.charAt(stop), sequence.charAt(stop+1))) {
				stop++;
			}
			for (int i=starts.size()-1; (i>starts.size()-1-totalAllowedStarts)&&i>=0; i--) {
				peptides.addAll(getPeptides(starts.get(i), stop, minLength, maxLength, sequence, constants));
			}
			starts.add(stop+1);
		}

		// digestion for sequence variants
		char stopCodon='*';
		int currentIndex=1;
		int stopCodonIndex;
		int cuts;
		int endIndex;
		int index;
		int start;
		TIntArrayList blockedIndices;
		HashMap<Integer, ArrayList<Integer>> usedPair=new HashMap<Integer, ArrayList<Integer>>();
		String sequenceVariant="";
		Collections.sort(variants);

		for (int i=0; i<variants.size(); i++) {
			AlleleVariant variant=variants.get(i);
			for (int idx=currentIndex; idx<starts.size(); idx++) {
				if ((variant.getStartSite()>starts.get(idx-1))&&(variant.getStartSite()<=starts.get(idx))) {
					currentIndex=idx;
					break;
				}
			}

			String perviousAA=(variant.getStartSite()-2<0)?"":sequence.substring(variant.getStartSite()-2, variant.getStartSite()-1);
			String nextAA=(variant.getStopSite()==sequence.length())?"*":sequence.substring(variant.getStopSite(), variant.getStopSite()+1);
			TIntArrayList addedStarts=getStartsAddedByVariant(variant, perviousAA, nextAA);

			blockedIndices=new TIntArrayList();
			if (variant.getStartSite()==starts.get(currentIndex)) {
				blockedIndices.add(currentIndex);
			} else if ((variant.getStartSite()-1==starts.get(currentIndex-1))&&!addedStarts.contains(starts.get(currentIndex-1))) {
				blockedIndices.add(currentIndex-1);
			}
			
			if (addedStarts.contains(starts.get(currentIndex-1))) {
				addedStarts.remove(starts.get(currentIndex-1));
			}
			
			
			endIndex=getNextIndex(currentIndex, 0, blockedIndices, starts.size()-1, 1);
			stopCodonIndex=variant.getNewSequence().indexOf(stopCodon);
			if (stopCodonIndex<0) {
				sequenceVariant=sequence.substring(0, variant.getStartSite()-1)+variant.getNewSequence()+sequence.substring(variant.getStopSite());
			} else {
				sequenceVariant=sequence.substring(0, variant.getStartSite()-1)+variant.getNewSequence().substring(0, stopCodonIndex);
			}

			usedPair=new HashMap<Integer, ArrayList<Integer>>();
			for (int j=0; j<addedStarts.size()+totalAllowedStarts; j++) {
				if (j<addedStarts.size()) {
					stop=addedStarts.get(j)-1;
				} else {
					index=getNextIndex(endIndex, j-addedStarts.size(), blockedIndices, starts.size()-1, 1);
					stop=starts.get(index)-1+variant.getNewSequence().length()-variant.getOriginalSequence().length();
				}
				if (stop>sequenceVariant.length()-1) {
					break;
				}

				cuts=totalAllowedStarts;
				while (cuts>0&&(j-cuts-addedStarts.size()<0)) {
					int offset=(j-cuts);
					if (offset<0) {
						index=getNextIndex(endIndex, 0-offset, blockedIndices, starts.size()-1, -1);
						start=starts.get(index);
					} else {
						start=addedStarts.get(offset);
					}
					if (!usedPair.containsKey(start)||!usedPair.get(start).contains(stop)) {
						peptides.addAll(getPeptides(start, stop, minLength, maxLength, sequenceVariant, constants));
						if (!usedPair.containsKey(start)) {
							usedPair.put(start, new ArrayList<Integer>());
						}
						usedPair.get(start).add(stop);
					}
					cuts--;
				}
			}
		}
		return peptides;
	}
	
	//@MoMo 
	private int getNextIndex(int index, int indexOffset, TIntArrayList blockedIndices, int lastIndex, int direction) {
		int nextIndex=index;
		while (indexOffset>=0) {
			nextIndex=(indexOffset!=0)?nextIndex+direction:nextIndex;
			while (blockedIndices.contains(nextIndex)) {
				nextIndex+=direction;
			}
			indexOffset--;
		}
		nextIndex=(nextIndex<0)?0:nextIndex;
		nextIndex=(nextIndex>lastIndex)?lastIndex:nextIndex;
		return nextIndex;
	}
	
	//@MoMo 
	private TIntArrayList getStartsAddedByVariant(AlleleVariant variant, String before, String after) {
		char stopCodon='*';
		TIntArrayList addedStarts=new TIntArrayList();
		String sequence=before+variant.getNewSequence()+after;
		int offset=before.length();
		for (int i=0; i<sequence.length()-1; i++) {
			if (isCutSite(sequence.charAt(i), sequence.charAt(i+1))) {
				addedStarts.add(i-offset+variant.getStartSite());
			} else if (sequence.charAt(i)==stopCodon) {
				addedStarts.add(i-offset+variant.getStartSite()-1);
				break;
			}
		}
		return addedStarts;
	}
	
	//@MoMo 
	private ArrayList<String> getPeptides(int start, int stop, int minLength, int maxLength, String sequence, AminoAcidConstants constants) {
		TCharDoubleHashMap fixedMods=constants.getFixedMods();
		ModificationMassMap variableMods=constants.getVariableMods();
		ArrayList<String> peptides=new ArrayList<String>();

		String peptide=sequence.substring(start, stop+1);
		if ((peptide.length()>=minLength)&&(peptide.length()<=maxLength)) {
			peptides.addAll(getModifiedForms(peptide, fixedMods, variableMods));

			if (start==0&&(variableMods!=null&&!variableMods.isEmpty()&&peptide.length()!=0)) {
				double mass=variableMods.getProteinNTermMod(peptide.charAt(0));
				if (mass!=ModificationMassMap.MISSING) {
					peptides.add("["+mass+"]"+peptide);
				}
			}
			if (stop==sequence.length()-1&&(variableMods!=null&&!variableMods.isEmpty()&&peptide.length()!=0)) {
				double mass=variableMods.getProteinCTermMod(peptide.charAt(peptide.length()-1));
				if (mass!=ModificationMassMap.MISSING) {
					peptides.add(peptide+"["+mass+"]");
				}
			}
		}
		return peptides;
	}



	public ArrayList<String> getModifiedForms(String peptide, TCharDoubleHashMap fixedMods, ModificationMassMap variableMods) {
		
		ArrayList<String> peptides=new ArrayList<String>();
		peptides.add(adjustForFixed(peptide, fixedMods));
		
		if (variableMods==null|| variableMods.isEmpty()||peptide.length()==0) return peptides;

		double mass=variableMods.getNTermMod(peptide.charAt(0));
		if (mass!=ModificationMassMap.MISSING) {
			peptides.add(adjustForFixed("["+mass+"]"+peptide, fixedMods));
		}
		mass=variableMods.getCTermMod(peptide.charAt(peptide.length()-1));
		if (mass!=ModificationMassMap.MISSING) {
			peptides.add(adjustForFixed(peptide+"["+mass+"]", fixedMods));
		}
		
		for (int i=0; i<peptide.length(); i++) {
			mass=variableMods.getVariableMod(peptide.charAt(i));
			if (mass!=ModificationMassMap.MISSING) {
				peptides.add(adjustForFixed(peptide.substring(0, i+1)+"["+mass+"]"+peptide.substring(i+1), fixedMods));
			}
		}
		return peptides;
	}
	public ArrayList<String> getTerminalModifiedForms(String peptide,int start,int end, int sequenceLength) {
		ArrayList<String> peptides=new ArrayList<String>();
		return peptides;
	}
	
	/**
	 * assumes if there's a mod in []s then the fixed part has already been considered
	 * @param peptide
	 * @param fixedMods
	 * @return
	 */
	public static String adjustForFixed(String peptide, TCharDoubleHashMap fixedMods) {
		StringBuilder sb=new StringBuilder();
		
		for (int i=0; i<peptide.length(); i++) {
			char aa=peptide.charAt(i);
			sb.append(aa);
			if (fixedMods.contains(aa)) {
				if (i==peptide.length()-1||peptide.charAt(i+1)!='[') {
					double mass=fixedMods.get(aa);
					sb.append('[');
					sb.append(mass);
					sb.append(']');
				}
			}
		}
		return sb.toString();
	}

}
