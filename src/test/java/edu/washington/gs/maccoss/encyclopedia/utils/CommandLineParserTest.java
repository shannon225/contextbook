package edu.washington.gs.maccoss.encyclopedia.utils;

import java.util.HashMap;
import java.util.Map.Entry;

import junit.framework.TestCase;

public class CommandLineParserTest extends TestCase {
	public static void main(String[] args) {
		HashMap<String, String> map=CommandLineParser.parseArguments(args);
		for (Entry<String, String> entry : map.entrySet()) {
			System.out.println(entry.getKey()+" = "+entry.getValue());
		}
		
	}
	public void testParsing() {
		String[] args=new String[] {"-q", "-i", "xxxxx", "-v"};
		HashMap<String, String> map=CommandLineParser.parseArguments(args);
		assertTrue(map.containsKey("-q"));
		assertTrue(map.containsKey("-i"));
		assertTrue(map.containsKey("-v"));
		assertEquals(null, map.get("-q"));
		assertEquals("xxxxx", map.get("-i"));
		assertEquals(null, map.get("-v"));
	}
}
