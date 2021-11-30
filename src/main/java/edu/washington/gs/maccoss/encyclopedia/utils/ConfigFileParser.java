package edu.washington.gs.maccoss.encyclopedia.utils;

import java.util.HashMap;
import java.util.Map.Entry;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class ConfigFileParser {
	public static final String CONFIG_FILE_TAG="-c";

	public static void updateArguments(HashMap<String, String> map) throws IOException {
		if (map.containsKey(CONFIG_FILE_TAG)) {
			InputStream input = new FileInputStream(map.get(CONFIG_FILE_TAG));
			Properties config = new Properties();
			config.load(input);
			
			for (Entry<Object, Object> e : config.entrySet()) {
				String key = (String) "-"+e.getKey();
				String value = (String) e.getValue();
				if (value.isEmpty()) {
					value = null;
				}
				if (!map.containsKey(key)) {
					map.put(key, value);
				}
			}
		}
	}

}
