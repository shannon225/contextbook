package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import org.apache.commons.lang3.ArrayUtils;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;

public class AmbiguousPeptideModSeq {
	private static final int NO_GROUP=0;
	
	private final String[] aas;
	private final boolean[] isModified;
	private final boolean[] modifiable;
	private final int[] modificationGroup;
	private final byte ambiguityDirection; // left is -1, right is 1, no direction is 0
	private final PeptideModification modification;
	
	private AmbiguousPeptideModSeq(String[] aas, boolean[] modifiable, boolean[] isModified, int[] modificationGroup, byte ambiguityDirection, PeptideModification modification) {
		this.aas=aas;
		this.modifiable=modifiable;
		this.isModified=isModified;
		this.modificationGroup=modificationGroup;
		this.ambiguityDirection=ambiguityDirection;
		this.modification=modification;
	}
	
	public int getAmbiguityValue() {
		int count=0;
		for (int i=0; i<modifiable.length; i++) {
			if (modifiable[i]&&modificationGroup[i]>0) count++;
		}
		return count;
	}
	
	public int getNumModifiableSites() {
		int count=0;
		for (int i=0; i<modifiable.length; i++) {
			if (modifiable[i]) count++;
		}
		return count;
	}
	
	public int getNumAmbigousResidues() {
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

	public Optional<AmbiguousPeptideModSeq> removeAmbiguity(PeptideModification modification, AmbiguousPeptideModSeq... confirmedIDs) {
		return removeAmbiguity(modification, Arrays.asList(confirmedIDs));
	}
	public Optional<AmbiguousPeptideModSeq> removeAmbiguity(PeptideModification modification, List<AmbiguousPeptideModSeq> confirmedIDs) {
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
		AmbiguousPeptideModSeq newSeq=new AmbiguousPeptideModSeq(newaas, newmodifiable, newisModified, newmodificationGroup, ambiguityDirection, modification);
		
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

	public String getNonDirectionalPeptideAnnotation() {
		return getPeptideAnnotation((byte)0);
	}
	public String getPeptideAnnotation() {
		return getPeptideAnnotation(ambiguityDirection);
	}
	public String getPeptideAnnotation(byte direction) {
		StringBuilder sb=new StringBuilder();
		
		int lastModGroup=-1;
		for (int i=0; i<aas.length; i++) {
			if (modifiable[i]) {
				if (modificationGroup[i]!=lastModGroup&&modificationGroup[i]!=0) {
					if (direction==-1) {
						sb.append("<");
					} else {
						sb.append("(");
					}
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
					if (direction==1) {
						sb.append(">");
					} else {
						sb.append(")");
					}
				}
			} else {
				sb.append(aas[i]);
			}
		}

		return sb.toString();
	}
	
	public static AmbiguousPeptideModSeq getAmbiguousPeptideModSeq(String peptideAnnotation, PeptideModification mod) {

		//String[] aas;
		//byte ambiguityDirection; // left is -1, right is 1, no direction is 0
		
		char[] ca=peptideAnnotation.toCharArray();
		
		ArrayList<String> aas=new ArrayList<String>();
		ArrayList<Boolean> isModified=new ArrayList<>();
		ArrayList<Boolean> modifiable=new ArrayList<>();
		TIntArrayList modificationGroup=new TIntArrayList();
		
		boolean inGroup=false;
		int groupNumber=0;
		byte ambiguityDirection=0;
		
		for (int i = 0; i < ca.length; i++) {
			if (ca[i]=='<'||ca[i]=='(') {
				// opening ambiguity
				groupNumber++;
				inGroup=true;
				ambiguityDirection=(ca[i]=='<')?(byte)-1:(byte)1;
				
			} else if (ca[i]=='>'||ca[i]==')') {
				// closing ambiguity
				inGroup=false;
				
			} else if (ca[i]=='[') {
				StringBuilder sb=new StringBuilder();
				i++;
				while (ca[i]!=']') {
					sb.append(ca[i]);
					i++;
				}
				if (aas.size()==0) {
					// pulls first amino acid in
					// this handling of n-termini mods assumes you can't have multiple []s in a row
					i++;
					aas.add(Character.toString(ca[i]));
				}
				String massText = sb.toString();
				double modificationMass = Double.valueOf(massText);
				String aaString=aas.get(aas.size()-1);
				char aaChar=aaString.charAt(0);
				modificationMass=MassConstants.getAccurateModificationMass(aaChar, modificationMass);

				// adjust last
				aas.set(aas.size()-1, aaString+(modificationMass>=0?"[+":"[")+modificationMass+"]");
				modificationGroup.removeAt(modificationGroup.size()-1);
				isModified.remove(isModified.size()-1);
				modifiable.remove(modifiable.size()-1);
				
				// add new
				isModified.add(mod.isModificationMass(aaChar, modificationMass));
				modifiable.add(mod.isModifiable(aaChar));
				if (inGroup) {
					modificationGroup.add(groupNumber);
				} else {
					modificationGroup.add(0);
				}
				
			} else {
				aas.add(Character.toString(ca[i]));
				isModified.add(false);
				modifiable.add(mod.isModifiable(ca[i]));
				if (inGroup) {
					modificationGroup.add(groupNumber);
				} else {
					modificationGroup.add(0);
				}
				
			}
		}
		
		String[] aaArray=aas.toArray(new String[aas.size()]);
		boolean[] isModifiedArray=ArrayUtils.toPrimitive(isModified.toArray(new Boolean[isModified.size()]));
		boolean[] isModifiableArray=ArrayUtils.toPrimitive(modifiable.toArray(new Boolean[modifiable.size()]));
		int[] modificationGroupArray=modificationGroup.toArray();
		
		return new AmbiguousPeptideModSeq(aaArray, isModifiableArray, isModifiedArray, modificationGroupArray, ambiguityDirection, mod);
	}
	
	public static AmbiguousPeptideModSeq getFullyAmbiguous(String targetPeptide, PeptideModification modification, AminoAcidConstants aaConstants) {
		String[] aas=PeptideUtils.getMasses(targetPeptide, aaConstants).z;
		boolean[] modifiable=new boolean[aas.length];
		boolean[] isModified=new boolean[aas.length];
		int[] modificationGroup=new int[aas.length];
		
		for (int i=0; i<aas.length; i++) {
			char c=aas[i].charAt(0);
			if (contains(modification.getModifiableAAs(), c)) {
				modifiable[i]=true;
				modificationGroup[i]=1;
			}
			if (modifiable[i]) {
				int mods=PeptideUtils.getNumberOfMods(aas[i], modification.getNominalMass());
				if (mods>0) {
					isModified[i]=true;
				}
			}
		}
		return new AmbiguousPeptideModSeq(aas, modifiable, isModified, modificationGroup, (byte)0, modification);
	}
	
	public static AmbiguousPeptideModSeq getUnambigous(String targetPeptide, PeptideModification modification, AminoAcidConstants aaConstants) {
		String[] aas=PeptideUtils.getMasses(targetPeptide, aaConstants).z;
		boolean[] modifiable=new boolean[aas.length];
		boolean[] isModified=new boolean[aas.length];
		int[] modificationGroup=new int[aas.length];
		
		int currentGroup=0;
		for (int i=0; i<aas.length; i++) {
			char c=aas[i].charAt(0);
			if (contains(modification.getModifiableAAs(), c)) {
				modifiable[i]=true;
			}
			if (modifiable[i]) {
				int mods=PeptideUtils.getNumberOfMods(aas[i], modification.getNominalMass());
				if (mods>0) {
					isModified[i]=true;
					currentGroup++;
					modificationGroup[i]=currentGroup;
				}
			}
		}
		return new AmbiguousPeptideModSeq(aas, modifiable, isModified, modificationGroup, (byte)0, modification);
	}
	
	public static AmbiguousPeptideModSeq getLeftAmbiguity(String targetPeptide, PeptideModification modification, AminoAcidConstants aaConstants) {
		String[] aas=PeptideUtils.getMasses(targetPeptide, aaConstants).z;
		boolean[] modifiable=new boolean[aas.length];
		boolean[] isModified=new boolean[aas.length];
		int[] modificationGroup=new int[aas.length];
		
		for (int i=0; i<aas.length; i++) {
			char c=aas[i].charAt(0);
			if (contains(modification.getModifiableAAs(), c)) {
				modifiable[i]=true;
			}
			if (modifiable[i]) {
				int mods=PeptideUtils.getNumberOfMods(aas[i], modification.getNominalMass());
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
		return new AmbiguousPeptideModSeq(aas, modifiable, isModified, modificationGroup, (byte)-1, modification);
	}
	
	public static AmbiguousPeptideModSeq getRightAmbiguity(String targetPeptide, PeptideModification modification, AminoAcidConstants aaConstants) {
		String[] aas=PeptideUtils.getMasses(targetPeptide, aaConstants).z;
		boolean[] modifiable=new boolean[aas.length];
		boolean[] isModified=new boolean[aas.length];
		int[] modificationGroup=new int[aas.length];
		
		for (int i=aas.length-1; i>=0; i--) {
			char c=aas[i].charAt(0);
			if (contains(modification.getModifiableAAs(), c)) {
				modifiable[i]=true;
			}
			if (modifiable[i]) {
				int mods=PeptideUtils.getNumberOfMods(aas[i], modification.getNominalMass());
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
		return new AmbiguousPeptideModSeq(aas, modifiable, isModified, modificationGroup, (byte)1, modification);
	}

	public static boolean isSiteSpecific(AmbiguousPeptideModSeq targetPeptideName, PeptideModification modification) {
		return isSiteSpecific(targetPeptideName.getPeptideAnnotation(), modification);
	}

	public static boolean isSiteSpecificAtEnd(AmbiguousPeptideModSeq targetPeptideName, PeptideModification modification) {
		return isSiteSpecificAtEnd(targetPeptideName.getPeptideAnnotation(), modification);
	}
	
	public static boolean isSiteSpecificAtEnd(String targetPeptideName, PeptideModification modification) {
		if (!isSiteSpecific(targetPeptideName, modification)) {
			return false;
		}

		char[] ca=targetPeptideName.toCharArray();
		
		int direction=0;
		for (int i = 0; i < ca.length; i++) {
			if (ca[i]=='<') {
				direction=-1;
				break;
			} else if (ca[i]=='>') {
				direction=1;
				break;
			}
		}
		
		int modsToFind=PeptideUtils.getNumberOfMods(targetPeptideName, modification.getNominalMass());
		if (direction==-1) {
			for (int i = 0; i < ca.length; i++) {
				if (contains(modification.getModifiableAAs(), ca[i])) {
					return false;
				}
				if (ca[i]=='<') {
					modsToFind--;
					while (i<ca.length&&ca[i]!=')') {
						i++;
					}
				}
				if (modsToFind==0) return true;
			}
			return true;
		}
		if (direction==1) {
			for (int i = ca.length-1; i >=0; i--) {
				if (contains(modification.getModifiableAAs(), ca[i])) {
					return false;
				}
				if (ca[i]=='>') {
					modsToFind--;
					while (i>=0&&ca[i]!='(') {
						i--;
					}
				}
				if (modsToFind==0) return true;
			}
			return true;
		}
		
		// can't determine direction!
		return false;
	}

	public static boolean isCompletelyAmbiguous(AmbiguousPeptideModSeq targetPeptideName, PeptideModification modification) {
		return isCompletelyAmbiguous(targetPeptideName.getPeptideAnnotation(), modification);
	}
	
	public static boolean isCompletelyAmbiguous(String targetPeptideName, PeptideModification modification) {
		char[] ca=targetPeptideName.toCharArray();

		StringBuilder sb=new StringBuilder();
		for (int i = 0; i < ca.length; i++) {
			if (ca[i]=='('||ca[i]=='<') {
				i++;
				while (ca[i]!=')'&&ca[i]!='>') {
					i++;
				}
			} else {
				sb.append(ca[i]);
			}
		}
		return getNumberOfSTYs(sb.toString(), modification.getModifiableAAs())==0;
	}
	
	public static boolean isSiteSpecific(String targetPeptideName, PeptideModification modification) {
		char[] ca=targetPeptideName.toCharArray();

		for (int i = 0; i < ca.length; i++) {
			if (ca[i]=='('||ca[i]=='<') {
				StringBuilder sb=new StringBuilder();
				i++;
				while (ca[i]!=')'&&ca[i]!='>') {
					sb.append(ca[i]);
					i++;
				}
				String massText = sb.toString();
				int mods=PeptideUtils.getNumberOfMods(massText, modification.getNominalMass());
				int modables=getNumberOfSTYs(massText, modification.getModifiableAAs());
				if (modables>mods) return false;
			}
		}
		return true;
	}
	
	private static int getNumberOfSTYs(String sequence, char[] modifiableAAs) {
		int total=0;
		for (int i=0; i<sequence.length(); i++) {
			char c=sequence.charAt(i);
			if (contains(modifiableAAs, c)) {
				total++;
			}
		}
		return total;
	}

	private static boolean contains(char[] aas, char aa) {
		for (int i=0; i<aas.length; i++) {
			if (aas[i]==aa) return true;
		}
		return false;
	}
}
