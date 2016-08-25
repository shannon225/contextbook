package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;

public class ProteinGroupQuantifier {
	private final HashMap<String, ProteinGroup> groups=new HashMap<String, ProteinGroup>();
	private final TObjectFloatHashMap<ProteinGroup> intensities=new TObjectFloatHashMap<ProteinGroup>();
	private final TObjectIntHashMap<ProteinGroup> numQuantPeptides=new TObjectIntHashMap<ProteinGroup>();

	public ProteinGroupQuantifier(ArrayList<ProteinGroup> groupList) {
		for (ProteinGroup group : groupList) {
			for (String accession : group.getEquivalentAccessions()) {
				groups.put(accession, group);
			}
		}
	}

	public boolean addIntensity(Collection<String> accessions, float intensity) {
		ProteinGroup group=null;
		for (String accession : accessions) {
			ProteinGroup newGroup=groups.get(accession);
			if (newGroup!=null) {
				if (group==null) {
					group=newGroup;
				} else if (group.equals(newGroup)) {
					continue; // still in this group
				} else {
					return false; // not unique!
				}
			}
		}
		if (group==null) return false;
		
		intensities.adjustOrPutValue(group, intensity, intensity);
		numQuantPeptides.adjustOrPutValue(group, 1, 1);
		return true;
	}
	
	public float getIntensity(ProteinGroup group) {
		return intensities.get(group);
	}
	
	public int getNumberOfQuantitativePeptides(ProteinGroup group) {
		return numQuantPeptides.get(group);
	}
}
