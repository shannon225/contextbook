package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;

public class CompoundFilenameFilter implements FilenameFilter {
	ArrayList<FilenameFilter> filters=new ArrayList<FilenameFilter>();

	public CompoundFilenameFilter(FilenameFilter... filter) {
		for (int i=0; i<filter.length; i++) {
			filters.add(filter[i]);
		}
	}
	
	public boolean accept(File dir, String name) {
		for (FilenameFilter filter : filters) {
			if (filter.accept(dir, name)) return true;
		}
		return false;
	};
}
