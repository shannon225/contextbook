package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import junit.framework.TestCase;

public class PercolatorReaderTest extends TestCase {
	private static final String REVERSE_PSMID="110415_bcs_hela_starved_DDA.mzML:11.096461:decoyEDIT[+79.966331]PEPR+2";
	private static final String FORWARD_PSMID="110415_bcs_hela_starved_DDA.mzML:11.096461:PEPT[+79.966331]IDER+2";
	
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), false);
	
	public static void main(String[] args) {
		File f=new File("/Volumes/BriansSSD/nick/110815_bcs_hela_6mz_600_700.mzML.encyclopedia.decoy.txt");
		ArrayList<PercolatorPeptide> passingPeptidesFromTSV=PercolatorReader.getPassingPeptidesFromTSV(f, 1.0f, true).x;
		System.out.println(passingPeptidesFromTSV.size());
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
		
		LibraryEntry entry=new LibraryEntry("", new HashSet<String>(), 518.73841, (byte)2, "PEPT[+79.966331]IDER", 1, 0.0f, 0.0f, massArray, intensityArray);
		LibraryEntry reverse=entry.getDecoy(PARAMETERS);

		String psmid=PercolatorPeptide.getPSMID(entry, 11.096461f, DUMMY_DIA_FILE);
		System.out.println(psmid);
		assertEquals(FORWARD_PSMID, psmid);

		String revpsmid=PercolatorPeptide.getPSMID(reverse, 11.096461f, DUMMY_DIA_FILE);
		System.out.println(revpsmid);
		assertEquals(REVERSE_PSMID, revpsmid);
	}
	
	public void testIsPSMIDDecoy() {
		assertFalse(PercolatorPeptide.isPSMIDDecoy(FORWARD_PSMID));
		assertTrue(PercolatorPeptide.isPSMIDDecoy(REVERSE_PSMID));
	}

	public void testGetPeptideSequence() {
		assertEquals("PEPT[+79.966331]IDER", PercolatorPeptide.getPeptideSequence(FORWARD_PSMID));
		assertEquals("EDIT[+79.966331]PEPR", PercolatorPeptide.getPeptideSequence(REVERSE_PSMID));
	}
	
	public void testGetCharge() {
		assertEquals((byte)2, PercolatorPeptide.getCharge(FORWARD_PSMID));
		assertEquals((byte)2, PercolatorPeptide.getCharge(REVERSE_PSMID));
		
	}
	public void testGetFile() {
		assertEquals("110415_bcs_hela_starved_DDA.mzML", PercolatorPeptide.getFile(FORWARD_PSMID));
		assertEquals("110415_bcs_hela_starved_DDA.mzML", PercolatorPeptide.getFile(REVERSE_PSMID));
	}

	private static final StripeFileInterface DUMMY_DIA_FILE = new StripeFileInterface() {
		@Override
		public HashMap<Range, Float> getRanges() {
			throw new UnsupportedOperationException();
		}

		@Override
		public void openFile(File userFile) throws IOException, SQLException {
			throw new UnsupportedOperationException();
		}

		@Override
		public ArrayList<PrecursorScan> getPrecursors(float minRT, float maxRT) throws IOException, SQLException, DataFormatException {
			throw new UnsupportedOperationException();
		}

		@Override
		public ArrayList<Stripe> getStripes(double targetMz, float minRT, float maxRT, boolean sqrt) throws IOException, SQLException {
			throw new UnsupportedOperationException();
		}

		@Override
		public ArrayList<Stripe> getStripes(Range targetMzRange, float minRT, float maxRT, boolean sqrt) throws IOException, SQLException {
			throw new UnsupportedOperationException();
		}

		@Override
		public float getTIC() throws IOException, SQLException {
			throw new UnsupportedOperationException();
		}

		@Override
		public float getGradientLength() throws IOException, SQLException {
			throw new UnsupportedOperationException();
		}

		@Override
		public void close() {
			throw new UnsupportedOperationException();
		}

		@Override
		public boolean isOpen() {
			throw new UnsupportedOperationException();
		}

		@Override
		public File getFile() {
			throw new UnsupportedOperationException();
		}

		@Override
		public String getOriginalFileName() {
			return "110415_bcs_hela_starved_DDA.mzML";
		}
	};
}
