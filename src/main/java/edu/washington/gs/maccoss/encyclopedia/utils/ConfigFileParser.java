package edu.washington.gs.maccoss.encyclopedia.utils;

import org.ini4j.Wini;
import org.ini4j.Profile.Section;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Set;
import java.util.Map.Entry;

public class ConfigFileParser {
	public static void updateArguments(HashMap<String, String> map) throws IOException {
		if (map.containsKey("-c")) {
			Wini ini = new Wini(new File(map.get("-c")));
			Set<Entry<String, Section>> sections = ini.entrySet();
			System.out.println("TEST");
			System.out.println(sections);

			for (Entry<String, Section> e : sections) {
				Section section = e.getValue();
			
				Set<Entry<String, String>> values = section.entrySet();
				for (Entry<String, String> e2 : values) {
					String key = "-"+e2.getKey();
					String value = e2.getValue();
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

}
