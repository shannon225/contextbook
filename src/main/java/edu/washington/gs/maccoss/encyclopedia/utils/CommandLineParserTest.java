package edu.washington.gs.maccoss.encyclopedia.utils;

import java.util.HashMap;

import junit.framework.TestCase;

public class CommandLineParserTest extends TestCase {
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
