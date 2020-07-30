package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import gnu.trove.map.hash.TCharDoubleHashMap;
import junit.framework.TestCase;

public class LibraryEntryModifierTest extends TestCase {
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(50), new MassTolerance(50),
			DigestionEnzyme.getEnzyme("trypsin"), false, true, false);

	public void testModifyModelAtEverySite() {
		String peptideModSeq="AQK[+42.0]HS[+79.96633]DSSLEER";
		FragmentationModel origModel=PeptideUtils.getPeptideModel(peptideModSeq, PARAMETERS.getAAConstants());
		
		TCharDoubleHashMap fixedMods=new TCharDoubleHashMap();
		fixedMods.put('K',8.014199);
		fixedMods.put('R',10.008269);
		
		// test adjustment
		FragmentationModel model=LibraryEntryModifier.modifyModelAtEverySite(peptideModSeq, fixedMods, false, PARAMETERS);
		
		assertEquals(origModel.getChargedMass((byte)1)+8.014199+10.008269, model.getChargedMass((byte)1), 0.000001);
		
		double[] masses=model.getMasses();
		double[] expected=new double[] {71.037114,128.058578,178.11972699999998,137.058912,166.998359,115.026943,87.032028,87.032028,113.084064,129.042593,129.042593,166.10938000000002};
		
		for (int i = 0; i < masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.00001);
		}
		
		// test replacement
		model=LibraryEntryModifier.modifyModelAtEverySite(peptideModSeq, fixedMods, true, PARAMETERS);
		
		assertEquals(origModel.getChargedMass((byte)1)+8.014199+10.008269-42.010565, model.getChargedMass((byte)1), 0.000001);
		
		masses=model.getMasses();
		expected=new double[] {71.037114,128.058578,136.109162,137.058912,166.998359,115.026943,87.032028,87.032028,113.084064,129.042593,129.042593,166.10938000000002};
		
		for (int i = 0; i < masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.00001);
		}
	}
	

	public void testModifyModelAtEachSite() {
		String peptideModSeq="AQK[+42.0]HS[+79.96633]DSSLEER";
		TCharDoubleHashMap fixedMods=new TCharDoubleHashMap();
		fixedMods.put('S',79.96633);

		ArrayList<FragmentationModel> models=LibraryEntryModifier.modifyModelAtEachSite(peptideModSeq, fixedMods, true, PARAMETERS);
		
		HashSet<String> expected=new HashSet<>();
		expected.add("AQK[+42.010565]HS[79.96633]DSSLEER");
		expected.add("AQK[+42.010565]HS[+79.966331]DS[79.96633]SLEER");
		expected.add("AQK[+42.010565]HS[+79.966331]DSS[79.96633]LEER");

		assertEquals(expected.size(), models.size());
		for (FragmentationModel eachModel : models) {
			assertTrue(expected.contains(eachModel.toString()));
			System.out.println(eachModel.toString()+" --> "+eachModel.getChargedMass((byte)1));
		}
	}
}
