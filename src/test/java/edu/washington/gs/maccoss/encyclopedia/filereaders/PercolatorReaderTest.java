package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import junit.framework.TestCase;

public class PercolatorReaderTest extends TestCase {
	public static void main(String[] args) {
		File f=new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_980_1140.mzML.pecan.txt.xml");
		ArrayList<ScoredObject<String>> data=PercolatorReader.getPassingPeptidesFromXML(f, 1.0f);
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

}
