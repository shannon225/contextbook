package edu.washington.gs.maccoss.encyclopedia.utils;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

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

	/**
	 * Parse the given command line arguments into a list of paths specified by {@code -i}
	 * and the remaining parameters wtih {@link #parseArguments(String[])}.
	 */
	public static Pair<List<Path>, HashMap<String, String>> parseFilesAndArguments(String[] args) {
		return parseFilesAndArguments(args, "-i");
	}

	/**
	 * Parse the given command line arguments into a list of paths specified by {@code filesFlag}
	 * and the remaining parameters wtih {@link #parseArguments(String[])}.
	 */
	public static Pair<List<Path>, HashMap<String, String>> parseFilesAndArguments(String[] args, String filesFlag) {
		final List<String> mutableArgs = Lists.newLinkedList(Arrays.asList(args));
		final List<Path> paths = Lists.newArrayListWithExpectedSize(1);

		for (Iterator<String> iter = mutableArgs.iterator(); iter.hasNext();) {
			if (Objects.equals(filesFlag, iter.next())) {
				iter.remove(); // pop the flag

				paths.add(Paths.get(iter.next()));

				iter.remove(); // pop the path
			}
		}

		return new Pair<>(
				ImmutableList.copyOf(paths),
				parseArguments(mutableArgs.toArray(new String[0]))
		);
	}
}
