package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;

import java.util.ArrayList;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import edu.washington.gs.maccoss.encyclopedia.algorithms.SearchTestSupport;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import junit.framework.TestCase;

public class PeptideEncodingTest extends TestCase {
	
	public void testINDArray() {
		// INCORRECT orientation for DL4J!
		float[] d=new float[] {1, 2, 3};
		INDArray array=Nd4j.create(d);
		assertEquals("[    1.0000,    2.0000,    3.0000]", array.toString());
		
		
		// CORRECT orientation for DL4J!
		float[][] td=new float[][] {new float[] {1f, 2f, 3f}};
		array=Nd4j.create(td);
		assertEquals("[[    1.0000,    2.0000,    3.0000]]", array.toString());
	}
	
	public void testEncoding() throws Exception {
		System.out.println("Java: "+System.getProperty("java.version"));
		LibraryFile library=SearchTestSupport.getResultLibrary();
		
    	SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
    	AminoAcidConstants aaConstants=parameters.getAAConstants();
    	ArrayList<LibraryEntry> entries=library.getAllEntries(false, aaConstants);
		
		for (LibraryEntry entry : entries) {
			PeptideEncoding encoding=new PeptideEncoding(entry, entry.getScanStartTime(), parameters);
			EncodedAminoAcid[] aas=EncodedAminoAcid.getAAs(entry.getPeptideModSeq(), parameters.getAAConstants());
			INDArray encodeInput = encoding.encodeInput(aas);
			INDArray encodeResult = encoding.encodeResult();
			
			System.out.println("Input: ");
			System.out.println(General.toString(encodeInput.toFloatVector()));

			System.out.println("\nOutput: ");
			System.out.println(encodeResult);
			
			break;
		}
	}

}
