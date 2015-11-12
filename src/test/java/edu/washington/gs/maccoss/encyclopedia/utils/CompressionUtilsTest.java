package edu.washington.gs.maccoss.encyclopedia.utils;

import java.io.IOException;
import java.util.zip.DataFormatException;
import java.util.zip.Deflater;

import junit.framework.TestCase;

public class CompressionUtilsTest extends TestCase {
	public void testCompression() throws IOException, DataFormatException {
		try {
		     // Encode a String into bytes
		     String inputString = "blahblahblahblahblahblahblahblahblahblahblahblah";
		     byte[] input = inputString.getBytes("UTF-8");

		     // Compress the bytes
		     byte[] output = new byte[100];
		     Deflater compresser = new Deflater();
		     compresser.setInput(input);
		     compresser.finish();
		     int compressedDataLength = compresser.deflate(output);
		     compresser.end();
		     
		     byte[] result=CompressionUtils.decompress(output, input.length);
		     String outputString = new String(result, "UTF-8");
		     
		     assertEquals(inputString, outputString);
		     assertTrue(compressedDataLength<input.length);
		     
		 } catch(java.io.UnsupportedEncodingException ex) {
		     // handle
		 } catch (java.util.zip.DataFormatException ex) {
		     // handle
		 }
	}
}
