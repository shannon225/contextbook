package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.util.Map;

public interface TableParserMuscle {
	public void processRow(Map<String, String> row);
	public void cleanup();
}
