package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class DIAFileReader implements StripeFileReaderInterface {
	
	@Override
	public boolean accept(File dir, String name) {
		return name.toLowerCase().endsWith(StripeFile.DIA_EXTENSION);
	}
	
	@Override
	public boolean canTryToReadFile(File f) {
		if (!f.exists()||!f.isFile()||!f.canRead()) return false; 
		return f.getName().toLowerCase().endsWith(StripeFile.DIA_EXTENSION);
	}

	@Override
	public StripeFileInterface readStripeFile(File f, SearchParameters parameters) {
		if (canTryToReadFile(f)) {
			return openDIAFile(f);
		} else {
			throw new EncyclopediaException("Can't read file type "+f.getAbsolutePath());
		}
	}

	public static StripeFileInterface openDIAFile(File f) {
		try {
			StripeFileInterface stripefile=new StripeFile();
			stripefile.openFile(f);
			return stripefile;
		} catch (IOException ioe) {
			Logger.errorLine("Unexpected exception reading DIA file: "+f.getName());
			Logger.errorException(ioe);
			throw new EncyclopediaException("Error reading DIA file!", ioe);
		} catch (SQLException sqle) {
			Logger.errorLine("Unexpected exception reading DIA file: "+f.getName());
			Logger.logException(sqle);
			throw new EncyclopediaException("Error reading DIA file!", sqle);
		}
	}
}
