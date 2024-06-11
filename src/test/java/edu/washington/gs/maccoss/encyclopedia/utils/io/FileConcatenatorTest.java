package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;

import junit.framework.TestCase;

public class FileConcatenatorTest extends TestCase {
	public void testConcat() throws Exception {
		File f1=writeTempFile("/retention_times_synthetic.txt");
		File f2=writeTempFile("/truncated.msp");
		File f3=writeTempFile("/pecan.feature.txt");
		File f4=writeTempFile("/nextprot2017_testPEFF1.0rc25_small.peff");
		
		File output = File.createTempFile("output_", ".test");
		output.deleteOnExit();
		
		FileConcatenator.saveConcatenatedFile(output, f1, f2, f3);
		FileConcatenator.saveConcatenatedFile(output, f4);

		BufferedReader in=new BufferedReader(new FileReader(output));
		String eachline;
		while ((eachline=in.readLine())!=null) {
			System.out.println(eachline);
		}
		in.close();
		
		int numFiles=FileConcatenator.getNumberOfSubFiles(output);
		assertEquals(4, numFiles);

		File extracted=FileConcatenator.extractFile(output, 1);
		assertEquals(f1.length(), extracted.length());
		extracted=FileConcatenator.extractFile(output, 2);
		assertEquals(f2.length(), extracted.length());
		extracted=FileConcatenator.extractFile(output, 3);
		assertEquals(f3.length(), extracted.length());
		extracted=FileConcatenator.extractFile(output, 4);
		assertEquals(f4.length(), extracted.length());
	}
	
	public static File writeTempFile(String fileResourceName) throws Exception {
		InputStream is=FileConcatenatorTest.class.getResourceAsStream(fileResourceName);
		return writeTempFile(is);
	}

	/**
	 * this forces the new lines to be consistent across platforms for the test
	 * @param is
	 * @return
	 * @throws Exception
	 */
	public static File writeTempFile(InputStream is) throws Exception {
		File output = File.createTempFile("output_", ".test");
		output.deleteOnExit();
		FileWriter fileStream=new FileWriter(output);
		BufferedWriter out=new BufferedWriter(fileStream);

		BufferedReader in=new BufferedReader(new InputStreamReader(is));

		String eachline;
		while ((eachline=in.readLine())!=null) {
			out.write(eachline);
			out.newLine();
		}
		in.close();
		
		out.close();
		
		return output;
	}
}
