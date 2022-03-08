package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.io.File;
import java.io.FilenameFilter;

public class SimpleFilenameFilter implements FilenameFilter {
	private final String[] ext;

	public SimpleFilenameFilter(String... extensions) {
		ext=new String[extensions.length];
		for (int i=0; i<extensions.length; i++) {
			ext[i]=extensions[i].toLowerCase();
		}
	}
	
	@Override
	public boolean accept(File dir, String name) {
		String lower=name.toLowerCase();
		for (String extension : ext) {
			if (lower.endsWith(extension)) {
				return true;
			}
		}
		return false;
	}

	public boolean accept(String name) {
		return accept(null, name);
	}

}
