package edu.washington.gs.maccoss.encyclopedia.utils;

import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;

import java.io.IOException;
import java.util.HashMap;
import java.nio.file.Files;
import java.nio.file.Path;

import org.apache.commons.io.FileUtils;
import org.junit.Before;
import junit.framework.TestCase;

public class ConfigFileParserTest extends TestCase {
	protected Path tmpDir;

	@Before
	public void setUp() throws Exception {
		tmpDir = Files.createTempDirectory(getName());
		FileUtils.forceDeleteOnExit(tmpDir.toFile());
	}

	public void testParsing() throws IOException {
		Path config_path = EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/config.properties", tmpDir, "config", ".properties");

		String[] args=new String[] {"-q", "xxxxx", "-c", config_path.toString(), "-v"};
		HashMap<String, String> map=CommandLineParser.parseArguments(args);
		assertTrue(map.containsKey("-q"));
		assertTrue(map.containsKey("-c"));
		assertTrue(map.containsKey("-v"));
		assertEquals("xxxxx", map.get("-q"));
		assertEquals(config_path.toString(), map.get("-c"));
		assertEquals(null, map.get("-v"));

		ConfigFileParser.updateArguments(map);
		assertTrue(map.containsKey("-q"));
		assertTrue(map.containsKey("-c"));
		assertTrue(map.containsKey("-v"));
		assertTrue(map.containsKey("-acquisition"));
		assertTrue(map.containsKey("-quiet"));
		assertTrue(map.containsKey("-expectedPeakWidth"));
		assertEquals("xxxxx", map.get("-q"));
		assertEquals(config_path.toString(), map.get("-c"));
		assertEquals(null, map.get("-v"));
		assertEquals("DIA", map.get("-acquisition"));
		assertEquals("25", map.get("-expectedPeakWidth"));
		assertEquals(null, map.get("-quiet"));
	}
}
