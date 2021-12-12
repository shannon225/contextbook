package edu.washington.gs.maccoss.encyclopedia.utils;

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

}
