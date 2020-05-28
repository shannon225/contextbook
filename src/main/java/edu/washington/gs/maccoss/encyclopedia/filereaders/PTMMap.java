package edu.washington.gs.maccoss.encyclopedia.filereaders;


import java.util.HashMap;
import java.util.List;

import uk.ac.ebi.pride.utilities.pridemod.ModReader;
import uk.ac.ebi.pride.utilities.pridemod.model.PTM;
import uk.ac.ebi.pride.utilities.pridemod.model.Specificity;
import uk.ac.ebi.pride.utilities.pridemod.model.Specificity.AminoAcid;
import uk.ac.ebi.pride.utilities.pridemod.model.Specificity.Position;

public class PTMMap {
	private static final ModReader modReader=ModReader.getInstance();
	private static final HashMap<String, PostTranslationalModification> cache=new HashMap<>();
	
	public static boolean exists(String accession) {
		return getPTM(accession)!=PostTranslationalModification.nothing;
	}

	public static String getName(String accession) {
		return getPTM(accession).getName();
	}

	public static PostTranslationalModification getPTM(String name, String specificity) {
		String accession=name+" ("+specificity+")";
		PostTranslationalModification ptm=cache.get(accession);
		if (ptm!=null) return ptm;
		
		synchronized (modReader) {
			// check again in case the cache has been populated by the time we get a lock
			ptm=cache.get(accession); 
			if (ptm!=null) return ptm;

			Position pos=Specificity.parsePositon(specificity);
			AminoAcid aa=AminoAcid.NONE;
			if (pos==Position.NONE) {
				aa=Specificity.parseAminoAcid(Character.toString(specificity.charAt(0)));
			}
			List<PTM> newPTMs=modReader.getPTMListByEqualName(name);

			//System.out.println(name+"|"+pos+"|"+aa+" --> "+newPTMs.size()+" possible");
			
			PTM found=null;
			ptmloop: for (PTM possible : newPTMs) {
				for (Specificity spec : possible.getSpecificityCollection()) {
					//System.out.println(spec.getPosition()+"=="+pos+", "+spec.getName()+"=="+aa);
					if (spec.getPosition()==pos&&spec.getName()==aa) {
						found=possible;
						break ptmloop;
					}
				}
			}
			if (found==null) {
				ptm=PostTranslationalModification.nothing;
				cache.put(accession, PostTranslationalModification.nothing);
			} else {
				ptm=new PostTranslationalModification(found);
				cache.put(accession, ptm);
			}
		}
		return ptm;
	}

	public static PostTranslationalModification getPTM(String accession) {
		PostTranslationalModification ptm=cache.get(accession);
		if (ptm!=null) return ptm;
		
		synchronized (modReader) {
			// check again in case the cache has been populated by the time we get a lock
			ptm=cache.get(accession); 
			if (ptm!=null) return ptm;
			
			PTM newPTM=modReader.getPTMbyAccession(accession);
			if (newPTM==null) {
				ptm=PostTranslationalModification.nothing;
				cache.put(accession, PostTranslationalModification.nothing);
			} else {
				ptm=new PostTranslationalModification(newPTM);
				cache.put(accession, ptm);
			}
		}
		return ptm;
	}
	
	public static double getDeltaMass(String accession) {
		return getPTM(accession).getDeltaMass();
	}
	
	public static class PostTranslationalModification {
		public static final PostTranslationalModification nothing=new PostTranslationalModification("nothing", 0.0);
		private final String name;
		private final double deltaMass;

		public PostTranslationalModification(PTM ptm) {
			name=ptm.getName();
			deltaMass=ptm.getMonoDeltaMass();
		}
		
		public PostTranslationalModification(String name, double deltaMass) {
			this.name=name;
			this.deltaMass=deltaMass;
		}

		public double getDeltaMass() {
			return deltaMass;
		}

		public String getName() {
			return name;
		}
	}
}
