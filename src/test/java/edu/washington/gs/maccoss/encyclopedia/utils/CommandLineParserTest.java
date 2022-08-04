package edu.washington.gs.maccoss.encyclopedia.utils;

import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;

import com.google.common.collect.ImmutableMap;
import junit.framework.TestCase;
import org.junit.Assert;

import static org.junit.Assert.assertArrayEquals;

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

	public void testUnparsing() {
		String[] expectedArgs=new String[] {"-q", "-i", "xxxxx", "-v"};
		HashMap<String, String> expectedMap=CommandLineParser.parseArguments(expectedArgs);

		HashMap<String, String> map=new HashMap<String, String>();
		map.put("-q", null);
		map.put("-i", "xxxxx");
		map.put("-v", null);
		String[] actualArgs = CommandLineParser.unparseArguments(map);
		HashMap<String, String> actualMap=CommandLineParser.parseArguments(actualArgs);
		
		assertTrue(expectedMap.equals(actualMap));
	}

	public void testParseMultiple() {
		String[] args = new String [] {"-a", "A", "-i", "one", "-b", "B", "-i", "two", "-c", "C"};

		final Pair<List<String>, List<String>> result = CommandLineParser.parseMultipleAndGetRemainingArguments(args, "-i");

		assertArrayEquals(
				new String[] {"one", "two"},
				result.x.toArray()
		);

		assertArrayEquals(
				new String[] {"-a", "A", "-b", "B", "-c", "C"},
				result.y.toArray()
		);
	}

	public void testParseMultipleAndRemaining() {
		String[] args = new String [] {"-a", "A", "-i", "one", "-b", "B", "-i", "two", "-c", "C"};

		final Pair<List<String>, HashMap<String, String>> result = CommandLineParser.parseMultipleAndRemainingArguments(args, "-i");

		assertArrayEquals(
				new String[] {"one", "two"},
				result.x.toArray()
		);

		assertEquals(
				ImmutableMap.of("-a", "A", "-b", "B", "-c", "C"),
				result.y
		);
	}
}
