package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class StripeFileGenerator {
	private static final StripeFileReaderInterface[] readers=new StripeFileReaderInterface[] {
		new MzmlToDIAConverter()	
	};

	public static StripeFileInterface getFile(File f, SearchParameters parameters) {
		if (!f.exists()||!f.canRead()) {
			throw new EncyclopediaException("Can't read file "+f.getAbsolutePath());
		}
		
		// first try to read if .DIA
		if (f.getName().toLowerCase().endsWith(StripeFile.DIA_EXTENSION)) {
			return openDIAFile(f);
		}
		
		// then try to change name to .DIA and read
		String absolutePath=f.getAbsolutePath();
		File diaFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+StripeFile.DIA_EXTENSION);
		if (diaFile.exists()&&diaFile.canRead()) {
			return openDIAFile(diaFile);
		}
		
		// otherwise try readers in order
		for (StripeFileReaderInterface reader : readers) {
			if (reader.canTryToReadFile(f)) {
				return reader.createStripeFile(f, parameters);
			}
		}
		
		throw new EncyclopediaException("Can't read file type "+f.getAbsolutePath());
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
