package edu.washington.gs.maccoss.encyclopedia.filereaders;

import junit.framework.TestCase;

public class PercolatorReaderTest extends TestCase {
	public void testParsing() {
		String eachline="<q_value>3.385e-01</q_value>";
		float f=Float.parseFloat(eachline.substring(9, eachline.length()-10));
		System.out.println(eachline.substring(9, eachline.length()-10)+" --> "+f);
		eachline="<q_value>0.000000e+00</q_value>";
		f=Float.parseFloat(eachline.substring(9, eachline.length()-10));
		System.out.println(eachline.substring(9, eachline.length()-10)+" --> "+f);
	}

}
