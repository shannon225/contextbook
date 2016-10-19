package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.set.hash.TIntHashSet;

public class AmbiguousPeptideModSeq {
	public static final int NOMINAL_MASS=80;
	private static final int NO_GROUP=0;
	
	private final String[] aas;
	private final boolean[] isModified;
	private final boolean[] modifiable;
	private final int[] modificationGroup;
	
	public AmbiguousPeptideModSeq(String[] aas, boolean[] modifiable, boolean[] isModified, int[] modificationGroup) {
		this.aas=aas;
		this.modifiable=modifiable;
		this.isModified=isModified;
		this.modificationGroup=modificationGroup;
	}
	
	public int numAmbigousResidues() {
		int currentIndex=0;
		int currentTotal=0;
		int currentGap=0;
		for (int i=0; i<modificationGroup.length; i++) {
			if (modificationGroup[i]>0) {
				if (modificationGroup[i]!=currentIndex) {
					currentIndex=modificationGroup[i];
				} else {
					currentTotal+=currentGap;
				}
				currentGap=0;
				currentTotal++;
				
			} else if (currentIndex>0) {
				currentGap++;
			}
		}
		return currentTotal;
	}
	
	public int length() {
		return aas.length;
	}
	
	int[] getModificationGroup() {
		return modificationGroup;
	}
	
	public TIntHashSet[] getAmbiguityGroups() {
		int groupCount=General.max(modificationGroup);
		TIntHashSet[] groups=new TIntHashSet[groupCount];
		for (int i=0; i<groups.length; i++) {
			groups[i]=new TIntHashSet();
		}
		for (int i=0; i<modificationGroup.length; i++) {
			if (modificationGroup[i]>NO_GROUP) {
				int index=modificationGroup[i]-1;
				groups[index].add(i);
			}
		}
		return groups;
	}
	
	public static int[] getModificationGroupsFromSets(TIntHashSet[] sets, int length) {
		int[] modificationGroup=new int[length];
		for (int i=0; i<sets.length; i++) {
			int setNumber=i+1;
			for (int index : sets[i].toArray()) {
				modificationGroup[index]=setNumber;
			}
		}
		return modificationGroup;
	}

	public Optional<AmbiguousPeptideModSeq> removeAmbiguity(AmbiguousPeptideModSeq... confirmedIDs) {
		return removeAmbiguity(Arrays.asList(confirmedIDs));
	}
	public Optional<AmbiguousPeptideModSeq> removeAmbiguity(List<AmbiguousPeptideModSeq> confirmedIDs) {
		String[] newaas=aas.clone();
		boolean[] newisModified=isModified.clone();
		boolean[] newmodifiable=modifiable.clone();
		int[] newmodificationGroup=modificationGroup.clone();
		
		for (AmbiguousPeptideModSeq ided : confirmedIDs) {
			TIntHashSet[] groups=ided.getAmbiguityGroups();
			if (groups.length>1) {
				// TODO think about how to tag IDed sites if there are more than one mod. This is tricky!
				// i.e. if you see:
				//            previous id: S(S[+80])QQQS(S[+80])R 
				//            in question: (SS[+80])QQQS(S[+80])R
				// then you know you have: (S[+80])SQQQS(S[+80])R
				//
				// but if you see:
				//            previous id: S(S[+80])QQQS(S[+80])R 
				//            in question: (SS[+80])QQQ(SS[+80])R
				// then you don't know anything new! The in question peptide could be: 
				//                 either: S(S[+80])QQQ(S[+80])SR 
				//                     or: (S[+80])SQQQS(S[+80])R
				//                     or: (S[+80])SQQQ(S[+80])SR
				return Optional.of(this);
			}
			if (groups[0].size()==1) {
				// TODO currently can only work with fully localized sites! Think about the complications here if we change this! (There are a lot of complications, this may be impossible)
				int[] values=groups[0].toArray();
				for (int i=0; i<values.length; i++) {
					int prevGroup=newmodificationGroup[values[i]];
					newmodificationGroup[values[i]]=NO_GROUP;
					newmodifiable[values[i]]=false;
					if (newisModified[values[i]]) {
						for (int j=0; j<newmodificationGroup.length; j++) {
							if (newmodificationGroup[j]==prevGroup) {
								// only move mod if we find another slot!
								newisModified[values[i]]=false;
								String mod=newaas[values[i]].substring(1);
								newaas[values[i]]=Character.toString(newaas[values[i]].charAt(0));
								newisModified[j]=true;
								newaas[j]=newaas[j]+mod;
								break;
							}
						}
					}
				}
			}
		}
		AmbiguousPeptideModSeq newSeq=new AmbiguousPeptideModSeq(newaas, newmodifiable, newisModified, newmodificationGroup);
		
		boolean ok=false;
		for (int i=0; i<newmodificationGroup.length; i++) {
			if (newmodificationGroup[i]>0) {
				ok=true;
				break;
			}
		}
		if (!ok) {
			// removed all ambiguity
			return Optional.ofNullable(null);
		}
		return Optional.of(newSeq);
	}
	
	@Override
	public int hashCode() {
		return getPeptideAnnotation().hashCode();
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj==null||!(obj instanceof AmbiguousPeptideModSeq)) {
			return false;
		}
		return getPeptideAnnotation().equals(((AmbiguousPeptideModSeq)obj).getPeptideAnnotation());
	}
	
	public String getPeptideModSeq() {
		return PeptideUtils.getSequence(aas);
	}
	
	public String toString() {
		StringBuilder sb1=new StringBuilder();
		StringBuilder sb2=new StringBuilder();
		StringBuilder sb3=new StringBuilder();
		StringBuilder sb4=new StringBuilder();
		for (int i=0; i<aas.length; i++) {
			sb1.append(General.formatCellToWidth(aas[i], 4));
			sb2.append(General.formatCellToWidth(""+isModified[i], 4));
			sb3.append(General.formatCellToWidth(""+modifiable[i], 4));
			sb4.append(General.formatCellToWidth(""+modificationGroup[i], 4));
		}
		return sb1.toString()+"\n"+sb2.toString()+"\n"+sb3.toString()+"\n"+sb4.toString();
	}
	
	public String getPeptideAnnotation() {
		StringBuilder sb=new StringBuilder();
		
		int lastModGroup=-1;
		for (int i=0; i<aas.length; i++) {
			if (modifiable[i]) {
				if (modificationGroup[i]!=lastModGroup&&modificationGroup[i]!=0) {
					sb.append("(");
					lastModGroup=modificationGroup[i];
				}
				sb.append(aas[i]);
				boolean endOfGroup=true;
				for (int j=i+1; j<aas.length; j++) {
					if (modificationGroup[j]==modificationGroup[i]) {
						endOfGroup=false;
					}
				}
				if (endOfGroup) {
					sb.append(")");
				}
			} else {
				sb.append(aas[i]);
			}
		}

		return sb.toString();
	}
	
	public static AmbiguousPeptideModSeq getFullyAmbiguous(String targetPeptide, AminoAcidConstants aaConstants) {
		String[] aas=PeptideUtils.getMasses(targetPeptide, aaConstants).z;
		boolean[] modifiable=new boolean[aas.length];
		boolean[] isModified=new boolean[aas.length];
		int[] modificationGroup=new int[aas.length];
		
		for (int i=0; i<aas.length; i++) {
			char c=aas[i].charAt(0);
			if (c=='S'||c=='T'||c=='Y') {
				modifiable[i]=true;
				modificationGroup[i]=1;
			}
			if (modifiable[i]) {
				int mods=PeptideUtils.getNumberOfMods(aas[i], NOMINAL_MASS);
				if (mods>0) {
					isModified[i]=true;
				}
			}
		}
		return new AmbiguousPeptideModSeq(aas, modifiable, isModified, modificationGroup);
	}
	
	public static AmbiguousPeptideModSeq getUnambigous(String targetPeptide, AminoAcidConstants aaConstants) {
		String[] aas=PeptideUtils.getMasses(targetPeptide, aaConstants).z;
		boolean[] modifiable=new boolean[aas.length];
		boolean[] isModified=new boolean[aas.length];
		int[] modificationGroup=new int[aas.length];
		
		int currentGroup=0;
		for (int i=0; i<aas.length; i++) {
			char c=aas[i].charAt(0);
			if (c=='S'||c=='T'||c=='Y') {
				modifiable[i]=true;
			}
			if (modifiable[i]) {
				int mods=PeptideUtils.getNumberOfMods(aas[i], NOMINAL_MASS);
				if (mods>0) {
					isModified[i]=true;
					currentGroup++;
					modificationGroup[i]=currentGroup;
				}
			}
		}
		return new AmbiguousPeptideModSeq(aas, modifiable, isModified, modificationGroup);
	}
	
	public static AmbiguousPeptideModSeq getLeftAmbiguity(String targetPeptide, AminoAcidConstants aaConstants) {
		String[] aas=PeptideUtils.getMasses(targetPeptide, aaConstants).z;
		boolean[] modifiable=new boolean[aas.length];
		boolean[] isModified=new boolean[aas.length];
		int[] modificationGroup=new int[aas.length];
		
		for (int i=0; i<aas.length; i++) {
			char c=aas[i].charAt(0);
			if (c=='S'||c=='T'||c=='Y') {
				modifiable[i]=true;
			}
			if (modifiable[i]) {
				int mods=PeptideUtils.getNumberOfMods(aas[i], NOMINAL_MASS);
				if (mods>0) {
					isModified[i]=true;
				}
			}
		}
		
		int currentGroup=1;
		for (int i=0; i<aas.length; i++) {
			if (modifiable[i]) {
				modificationGroup[i]=currentGroup;
				if (isModified[i]) {
					currentGroup++;
				}
				boolean keepAnnotating=false;
				for (int j=i+1; j<isModified.length; j++) {
					if (isModified[j]) {
						keepAnnotating=true;
						break;
					}
				}
				if (!keepAnnotating) {
					break;
				}
			}
		}
		return new AmbiguousPeptideModSeq(aas, modifiable, isModified, modificationGroup);
	}
	
	public static AmbiguousPeptideModSeq getRightAmbiguity(String targetPeptide, AminoAcidConstants aaConstants) {
		String[] aas=PeptideUtils.getMasses(targetPeptide, aaConstants).z;
		boolean[] modifiable=new boolean[aas.length];
		boolean[] isModified=new boolean[aas.length];
		int[] modificationGroup=new int[aas.length];
		
		for (int i=aas.length-1; i>=0; i--) {
			char c=aas[i].charAt(0);
			if (c=='S'||c=='T'||c=='Y') {
				modifiable[i]=true;
			}
			if (modifiable[i]) {
				int mods=PeptideUtils.getNumberOfMods(aas[i], NOMINAL_MASS);
				if (mods>0) {
					isModified[i]=true;
				}
			}
		}
		
		int currentGroup=1;
		for (int i=aas.length-1; i>=0; i--) {
			if (modifiable[i]) {
				modificationGroup[i]=currentGroup;
				if (isModified[i]) {
					currentGroup++;
				}
				boolean keepAnnotating=false;
				for (int j=i-1; j>=0; j--) {
					if (isModified[j]) {
						keepAnnotating=true;
						break;
					}
				}
				if (!keepAnnotating) {
					break;
				}
			}
		}
		return new AmbiguousPeptideModSeq(aas, modifiable, isModified, modificationGroup);
	}

	public static boolean isLocalized(AmbiguousPeptideModSeq targetPeptideName) {
		return isLocalized(targetPeptideName.getPeptideAnnotation());
	}
	public static boolean isLocalized(String targetPeptideName) {
		char[] ca=targetPeptideName.toCharArray();

		for (int i = 0; i < ca.length; i++) {
			if (ca[i]=='(') {
				StringBuilder sb=new StringBuilder();
				i++;
				while (ca[i]!=')') {
					sb.append(ca[i]);
					i++;
				}
				String massText = sb.toString();
				int mods=PeptideUtils.getNumberOfMods(massText, NOMINAL_MASS);
				int modables=getNumberOfSTYs(massText);
				if (modables>mods) return false;
			}
		}
		return true;
	}
	
	private static int getNumberOfSTYs(String sequence) {
		int total=0;
		for (int i=0; i<sequence.length(); i++) {
			char c=sequence.charAt(i);
			if (c=='S'||c=='T'||c=='Y') {
				total++;
			}
		}
		return total;
	}

}
