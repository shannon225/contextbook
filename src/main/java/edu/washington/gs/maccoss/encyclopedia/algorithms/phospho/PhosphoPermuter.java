package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import gnu.trove.list.array.TIntArrayList;

public class PhosphoPermuter {
	
	public static ArrayList<String> getPermutations(String peptideModSeq, PeptideModification modification, AminoAcidConstants aaConstants) {
		String modificationTag=modification.toMassString();
		char[] modifibleAAs=modification.getModifiableAAs();
		
		FragmentationModel model=PeptideUtils.getPeptideModel(peptideModSeq, aaConstants);
		String[] aas=model.getAas();
		
		int phosphoCount=0;
		TIntArrayList styIndices=new TIntArrayList();
		StringBuilder sb=new StringBuilder();
		for (int i=0; i<aas.length; i++) {
			Pair<Character, Double> aa=FragmentationModel.parseAA(aas[i]);
			char aaValue=aa.x.charValue();
			
			boolean isSTY=false;
			for (int j=0; j<modifibleAAs.length; j++) {
				if (modifibleAAs[j]==aaValue) {
					isSTY=true;
				}
			}
			
			if (isSTY) {
				if (aa.y!=null) {
					if (modification.isModificationMass(aaValue, aa.y)) {
						phosphoCount++;
					} else {
						// STY with a non-phospho mod, so add as if it were another amino acid and continue
						sb.append(aas[i]);
						continue;
					}
				}
				styIndices.add(sb.length());
				sb.append(aaValue);
			} else {
				// non-phospho mods get added as is
				sb.append(aas[i]);
			}
		}
		
		ArrayList<String> sequences=new ArrayList<String>();
		if (phosphoCount==0) {
			StringBuilder seq=new StringBuilder(sb);
			sequences.add(seq.toString());
			
		} else if (phosphoCount==1) {
			for (int j=0; j<styIndices.size(); j++) {
				StringBuilder seq=new StringBuilder(sb);
				seq.insert(styIndices.get(j)+1, modificationTag);
				sequences.add(seq.toString());
			}

		} else if (phosphoCount==2) {
			for (int j=0; j<styIndices.size(); j++) {
				for (int k=j+1; k<styIndices.size(); k++) {
					StringBuilder seq=new StringBuilder(sb);
					seq.insert(styIndices.get(j)+1, modificationTag);
					seq.insert(styIndices.get(k)+1+modificationTag.length(), modificationTag);
					sequences.add(seq.toString());
				}
			}
		} else if (phosphoCount==3) {
			for (int j=0; j<styIndices.size(); j++) {
				for (int k=j+1; k<styIndices.size(); k++) {
					for (int l=k+1; l<styIndices.size(); l++) {
						StringBuilder seq=new StringBuilder(sb);
						seq.insert(styIndices.get(j)+1, modificationTag);
						seq.insert(styIndices.get(k)+1+modificationTag.length(), modificationTag);
						seq.insert(styIndices.get(l)+1+modificationTag.length()+modificationTag.length(), modificationTag);
						sequences.add(seq.toString());
					}
				}
			}
		} else if (phosphoCount==4) {
			for (int j=0; j<styIndices.size(); j++) {
				for (int k=j+1; k<styIndices.size(); k++) {
					for (int l=k+1; l<styIndices.size(); l++) {
						for (int m=l+1; m<styIndices.size(); m++) {
							StringBuilder seq=new StringBuilder(sb);
							seq.insert(styIndices.get(j)+1, modificationTag);
							seq.insert(styIndices.get(k)+1+modificationTag.length(), modificationTag);
							seq.insert(styIndices.get(l)+1+modificationTag.length()+modificationTag.length(), modificationTag);
							seq.insert(styIndices.get(m)+1+modificationTag.length()+modificationTag.length()+modificationTag.length(), modificationTag);
							sequences.add(seq.toString());
						}
					}
				}
			}
		} else if (phosphoCount==5) {
			for (int j=0; j<styIndices.size(); j++) {
				for (int k=j+1; k<styIndices.size(); k++) {
					for (int l=k+1; l<styIndices.size(); l++) {
						for (int m=l+1; m<styIndices.size(); m++) {
							for (int n=m+1; n<styIndices.size(); n++) {
								StringBuilder seq=new StringBuilder(sb);
								seq.insert(styIndices.get(j)+1, modificationTag);
								seq.insert(styIndices.get(k)+1+modificationTag.length(), modificationTag);
								seq.insert(styIndices.get(l)+1+modificationTag.length()+modificationTag.length(), modificationTag);
								seq.insert(styIndices.get(m)+1+modificationTag.length()+modificationTag.length()+modificationTag.length(), modificationTag);
								seq.insert(styIndices.get(n)+1+modificationTag.length()+modificationTag.length()+modificationTag.length()+modificationTag.length(), modificationTag);
								sequences.add(seq.toString());
							}
						}
					}
				}
			}
			
		} else {
			throw new EncyclopediaException("Sorry, no support for more than 5 phosphorylations: (you tried "+phosphoCount+" from "+peptideModSeq+"). Please report back so we know how to deal with this in the future!");
		}
		return sequences;
	}
}
