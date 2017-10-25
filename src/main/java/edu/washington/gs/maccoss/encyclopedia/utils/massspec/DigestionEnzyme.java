package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.AlleleVariant;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.ExtendedFastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntry;
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
		String originalSeq=entry.getSequence();
		ArrayList<String> peptides=digestProtein(originalSeq, minLength, maxLength, maxMissedCleavages, constants);
		if (entry instanceof ExtendedFastaEntry) {
			ExtendedFastaEntry peffentry=(ExtendedFastaEntry)entry;
			ArrayList<AlleleVariant> variants=peffentry.getPotentialVariant();
			for (int index=0; index<variants.size(); index++) {
				AlleleVariant variant=variants.get(index);
				String newSeq=originalSeq.substring(0, variant.getStartSite()-1)+variant.getNewSequence()+originalSeq.substring(variant.getStopSite());
				peptides.addAll(this.digestProtein(newSeq, minLength, maxLength, maxMissedCleavages, constants));
			}
			peptides=new ArrayList<String>(new HashSet<String>(peptides));
		}
		return peptides;
	}
	
	private ArrayList<String> digestProtein(String sequence, int minLength, int maxLength, int maxMissedCleavages, AminoAcidConstants constants) {
		TCharDoubleHashMap fixedMods=constants.getFixedMods();
		ModificationMassMap variableMods=constants.getVariableMods();
		int totalAllowedStarts=maxMissedCleavages+1;
		ArrayList<String> peptides=new ArrayList<String>();
		String peptide;
		char stopCodon = '*';
		TIntArrayList starts=new TIntArrayList();
		starts.add(0);
		int stop;
		
		//@MoMo; adding stop codon checking
		while (starts.get(0)<sequence.length()&&sequence.charAt(starts.get(0))!=stopCodon) {
			stop=starts.get(0);
			while ((stop<sequence.length()-1)&&!isCutSite(sequence.charAt(stop), sequence.charAt(stop+1))&&sequence.charAt(stop+1)!=stopCodon) {
				stop++;
			}
			for (int i=0; i<starts.size(); i++) {
				int start=starts.get(i);
				peptide=sequence.substring(start, stop+1);
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
			}
			starts.insert(0, stop+1);
			if (starts.size()>totalAllowedStarts) {
				starts.removeAt(starts.size()-1);
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
