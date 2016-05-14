package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ReverseLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import junit.framework.TestCase;

public class PercolatorReaderTest extends TestCase {
	private static final String REVERSE_PSMID="110415_bcs_hela_starved_DDA.mzML:11.096461:decoyEDIT[+80.0]PEPR+2";
	private static final String FORWARD_PSMID="110415_bcs_hela_starved_DDA.mzML:11.096461:PEPT[+80]IDER+2";
	
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"));
	
	public static void main(String[] args) {
		File f=new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110315_bcs_hela_phospho_starved_DDA.dia.percolator.txt");
		ArrayList<ScoredObject<String>> data=PercolatorReader.getPassingPeptidesFromTSV(f, 1.0f);
		System.out.println(data.size());
	}
	public void testParsing() {
		String eachline="<q_value>3.385e-01</q_value>";
		float f=Float.parseFloat(eachline.substring(9, eachline.length()-10));
		System.out.println(eachline.substring(9, eachline.length()-10)+" --> "+f);
		eachline="<q_value>0.000000e+00</q_value>";
		f=Float.parseFloat(eachline.substring(9, eachline.length()-10));
		System.out.println(eachline.substring(9, eachline.length()-10)+" --> "+f);
	}

	public void testGetPSMID() {
		double[] massArray = new double[] { 98.06063, 175.11955, 227.10323,
				304.16214, 324.15599, 333, 419.18908, 444, 505.20367,
				532.27314, 555, 618.28773, 650.407259, 666, 713.32082,
				733.31467, 777, 779.449849, 810.37359, 862.35727, 876.502609,
				888, 939.41618, 1018.45838, 1036.46894 };
		float[] intensityArray = new float[] { 1f, 2f, 3f, 4f, 5f, 6f,
				7f, 8f, 9f, 10f, 11f, 12f, 13f, 14f, 15f, 16f, 17f, 18f, 19f,
				20f, 21f, 22f, 23f, 24f, 25f };
		
		LibraryEntry entry=new LibraryEntry("", new HashSet<String>(), 518.73841, (byte)2, "PEPT[+80]IDER", 1, 0.0f, 0.0f, massArray, intensityArray);
		ReverseLibraryEntry reverse=entry.getDecoy(PARAMETERS, false);

		File diaFile=new File("/Users/searleb/Documents/freezer_experiment/110815_hela_experiment/data/hela_experiment/110415_bcs_hela_starved_DDA.mzML"); // FIXME unit test is not platform independent (will fail on windows machines)
		String psmid=PercolatorReader.getPSMID(entry, 11.096461f, diaFile);
		System.out.println(psmid);
		assertEquals(FORWARD_PSMID, psmid);

		String revpsmid=PercolatorReader.getPSMID(reverse, 11.096461f, diaFile);
		System.out.println(revpsmid);
		assertEquals(REVERSE_PSMID, revpsmid);
	}
	
	public void testIsPSMIDDecoy() {
		assertFalse(PercolatorReader.isPSMIDDecoy(FORWARD_PSMID));
		assertTrue(PercolatorReader.isPSMIDDecoy(REVERSE_PSMID));
	}

	public void testGetPeptideSequence() {
		assertEquals("PEPT[+80]IDER", PercolatorReader.getPeptideSequence(FORWARD_PSMID));
		assertEquals("EDIT[+80.0]PEPR", PercolatorReader.getPeptideSequence(REVERSE_PSMID));
	}
	
	public void testGetCharge() {
		assertEquals((byte)2, PercolatorReader.getCharge(FORWARD_PSMID));
		assertEquals((byte)2, PercolatorReader.getCharge(REVERSE_PSMID));
		
	}
	public void testGetFile() {
		assertEquals("110415_bcs_hela_starved_DDA.mzML", PercolatorReader.getFile(FORWARD_PSMID));
		assertEquals("110415_bcs_hela_starved_DDA.mzML", PercolatorReader.getFile(REVERSE_PSMID));
	}
}
