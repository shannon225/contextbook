package edu.washington.gs.maccoss.encyclopedia.datastructures;

import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;

public class ProteinGroupQuantifier {
	private final HashMap<String, ProteinGroupInterface> groups=new HashMap<String, ProteinGroupInterface>();
	private final TObjectFloatHashMap<ProteinGroupInterface> intensities=new TObjectFloatHashMap<ProteinGroupInterface>();
	private final TObjectIntHashMap<ProteinGroupInterface> numQuantPeptides=new TObjectIntHashMap<ProteinGroupInterface>();

	public ProteinGroupQuantifier(List<? extends ProteinGroupInterface> groupList) {
		for (ProteinGroupInterface group : groupList) {
			for (String accession : group.getEquivalentAccessions()) {
				groups.put(accession, group);
			}
		}
	}

	public boolean addIntensity(Collection<String> accessions, float intensity) {
		ProteinGroupInterface group=null;
		for (String accession : accessions) {
			ProteinGroupInterface newGroup=groups.get(accession);
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
	
	public float getIntensity(ProteinGroupInterface group) {
		return intensities.get(group);
	}
	
	public int getNumberOfQuantitativePeptides(ProteinGroupInterface group) {
		return numQuantPeptides.get(group);
	}
}
