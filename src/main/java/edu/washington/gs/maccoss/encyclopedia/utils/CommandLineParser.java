package edu.washington.gs.maccoss.encyclopedia.utils;

import java.util.ArrayList;
import java.util.HashMap;

public class CommandLineParser {
	
	public static HashMap<String, String> parseArguments(String[] args) {
		HashMap<String, String> map=new HashMap<String, String>();
		if (args.length==0) {
			return map;
		}
		
		for (int i=0; i<args.length; i++) {
			if (i<args.length-1&&!args[i+1].startsWith("-")) {
				map.put(args[i], args[i+1]);
				i++;
			} else {
				map.put(args[i], null);
			}
		}
		return map;
	}

	public static String[] unparseArguments(HashMap<String, String> map) {
		ArrayList<String> args = new ArrayList<String>();
		if (map.isEmpty()) {
			return args.toArray(new String[0]);
		}
		
		for (HashMap.Entry<String, String> entry : map.entrySet()) {
			String key = entry.getKey();
			String value = entry.getValue();
			args.add(key);
			if (value != null) {
				args.add(value);
			}
		}
		return args.toArray(new String[0]);
	}

}
