package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.map.hash.TFloatObjectHashMap;
import gnu.trove.procedure.TFloatObjectProcedure;
import junit.framework.TestCase;

public class AminoAcidConstantsTest extends TestCase {
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(50), new MassTolerance(50),
			DigestionEnzyme.getEnzyme("trypsin"), false, true, false);
	
	public void testGetMass() {
		String pep="AAGPLLTDEC[+57.0]R";
		System.out.println(PARAMETERS.getAAConstants().getChargedMass(pep, (byte)2));
		System.out.println(PARAMETERS.getAAConstants().getMass('C'));
		System.out.println(PARAMETERS.getAAConstants().getMass(pep));
	}
	
	public void testPeptideModSeq() {
		String pep="AAGPLLTDECR";
		assertEquals("AAGPLLTDEC[+57.0214635]R", PARAMETERS.getAAConstants().toPeptideModSeq(pep));
	}
	
	public void testGetChargedMass() {
		AminoAcidConstants NO_MODS = new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());
		String pep="MAREC[+57.0214635]YSPEPTIDESK";
		assertEquals(956.9244, NO_MODS.getChargedMass(pep, (byte)2), 0.0001);
		
		AminoAcidConstants C_MODS = new AminoAcidConstants();
		assertEquals(956.9244, C_MODS.getChargedMass(pep, (byte)2), 0.0001);
	}
	
	public static void main(String[] args) {
		//delta masses
		AminoAcidConstants aas = new AminoAcidConstants();
		char[] AAs="ALVIHYWFMSTNDQEKRCGP".toCharArray();

		for (int i = 0; i < AAs.length; i++) {
			char a = AAs[i];
			float massA = Math.round(aas.getMass(a)*100f)/100f;

			for (int j = 0; j < AAs.length; j++) {
				if (j>0) System.out.print('\t');
				char b = AAs[j];
				float mass=Math.round(Math.abs(aas.getMass(a)-aas.getMass(b))*100)/100f;
				System.out.print(mass);
			}
			System.out.println();
		}
		
	}
	
	public static void main2(String[] args) {
		AminoAcidConstants aas = new AminoAcidConstants();
		char[] AAs="ARNDCEQGHLKMFPSTWYV".toCharArray();
		TFloatObjectHashMap<ArrayList<String>> massMap=new TFloatObjectHashMap<>();
		
		for (int i = 0; i < AAs.length; i++) {
			char a = AAs[i];
			float massA = Math.round(aas.getMass(a)*100f)/100f;
			ArrayList<String> list=massMap.get(massA);
			if (list==null) {
				list=new ArrayList<String>();
				massMap.put(massA, list);
			}
			list.add(Character.toString(a));
			

			for (int j = i; j < AAs.length; j++) {
				char b = AAs[j];
				float mass=Math.round((aas.getMass(a)+aas.getMass(b))*100)/100f;

				list=massMap.get(mass);
				if (list==null) {
					list=new ArrayList<String>();
					massMap.put(mass, list);
				}
				list.add(Character.toString(a)+Character.toString(b));
			}
		}
		massMap.forEachEntry(new TFloatObjectProcedure<ArrayList<String>>() {
			@Override
			public boolean execute(float a, ArrayList<String> b) {
				if (b.size()>1) System.out.println(a);
				return true;
			}
		});

//		System.out.print("#\t-");
//		for (char a : AAs) {
//			System.out.print("\t"+a);
//		}
//		System.out.println();
//
//		for (char a : AAs) {
//			System.out.print(a);
//			float mass=Math.round((aas.getMass(a))*100)/100f;
//			System.out.print("\t"+(massMap.get(mass).size()>1));
//			//System.out.print("\t"+mass);
//			
//			for (char b : AAs) {
//				mass=Math.round((aas.getMass(a)+aas.getMass(b))*100)/100f;
//				System.out.print("\t"+(massMap.get(mass).size()>1));
//				//System.out.print("\t"+mass);
//			}
//			System.out.println();
//		}
	}
}
