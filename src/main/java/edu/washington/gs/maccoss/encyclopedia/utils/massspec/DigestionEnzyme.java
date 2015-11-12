package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TCharHashSet;

import java.util.ArrayList;

public class DigestionEnzyme {
	private static final char[] AAs="ACDEFGHIKLMNPQRSTVWY".toCharArray();
	private final String name;
	private final TCharHashSet nterm;
	private final TCharHashSet cterm;
	
	public static DigestionEnzyme getEnzyme(String enzymeName) {
		TCharHashSet n=new TCharHashSet();
		TCharHashSet c=new TCharHashSet();
		if ("Trypsin".equalsIgnoreCase(enzymeName)) {
			n.add('K');
			n.add('R');
			c.addAll(AAs);
			c.remove('P');
			
			return new DigestionEnzyme("Trypsin", n, c);
			
		} else if ("Trypsin/p".equalsIgnoreCase(enzymeName)) {
			n.add('K');
			n.add('R');
			c.addAll(AAs);
			
			return new DigestionEnzyme("Trypsin/p", n, c);
			
		} else if ("Lys-C".equalsIgnoreCase(enzymeName)) {
			n.add('K');
			c.addAll(AAs);
			c.remove('P');
			
			return new DigestionEnzyme("Lys-C", n, c);
			
		} else if ("Lys-N".equalsIgnoreCase(enzymeName)) {
			n.addAll(AAs);
			c.add('K');
			
			return new DigestionEnzyme("Lys-N", n, c);
			
		} else if ("Arg-C".equalsIgnoreCase(enzymeName)) {
			n.add('R');
			c.addAll(AAs);
			c.remove('P');
			
			return new DigestionEnzyme("Arg-C", n, c);
			
		} else if ("CNBr".equalsIgnoreCase(enzymeName)) {
			n.add('M');
			c.addAll(AAs);
			
			return new DigestionEnzyme("CNBr", n, c);
			
		} else if ("Chymotrypsin".equalsIgnoreCase(enzymeName)) {
			n.add('F');
			n.add('Y');
			n.add('W');
			n.add('L');
			c.addAll(AAs);
			c.remove('P');
			
			return new DigestionEnzyme("Chymotrypsin", n, c);
			
		} else if ("PepsinA".equalsIgnoreCase(enzymeName)) {
			n.add('F');
			n.add('L');
			c.addAll(AAs);
			
			return new DigestionEnzyme("PepsinA", n, c);
		}
		
		throw new IllegalArgumentException("Unknown digestion enzyme ["+enzymeName+"]");
	}
	
	DigestionEnzyme(String name, TCharHashSet nterm, TCharHashSet cterm) {
		this.name=name;
		this.nterm=nterm;
		this.cterm=cterm;
	}
	
	public String getName() {
		return name;
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
	
	public ArrayList<String> digestProtein(String sequence, int minLength, int maxLength, int maxMissedCleavages) {
		ArrayList<String> peptides=new ArrayList<String>();
		String peptide;
		TIntArrayList starts=new TIntArrayList();
		starts.add(0);
		int stop;

		while (starts.get(0)<sequence.length()) {
			stop=starts.get(0);
			while ((stop<sequence.length()-1)&&!isCutSite(sequence.charAt(stop), sequence.charAt(stop+1))) {
				stop++;
			}
			for (int i=0; i<starts.size(); i++) {
				int start=starts.get(i);
				peptide=sequence.substring(start, stop+1);
				if ((peptide.length()>=minLength)&&(peptide.length()<=maxLength)) {
					peptides.add(peptide);
				}
			}
			starts.insert(0, stop+1);
			if (starts.size()>maxMissedCleavages) {
				starts.removeAt(starts.size()-1);
			}
		}
		return peptides;
	}

}
