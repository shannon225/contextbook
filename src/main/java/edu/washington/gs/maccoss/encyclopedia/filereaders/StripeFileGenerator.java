package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.FilenameFilter;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.gui.general.CompoundFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class StripeFileGenerator {
	private static final DIAFileReader DIA_FILE_READER=new DIAFileReader();
	private static final StripeFileReaderInterface[] readers=new StripeFileReaderInterface[] {
			DIA_FILE_READER,
			new MzmlToDIAConverter()
	};

	public static StripeFileInterface getFile(File f, SearchParameters parameters) {
		// try to change name to .DIA and read
		String absolutePath=f.getAbsolutePath();
		File diaFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+StripeFile.DIA_EXTENSION);
		if (diaFile.exists()&&diaFile.canRead()) {
			try {
				return DIA_FILE_READER.readStripeFile(diaFile, parameters);
			} catch (EncyclopediaException ee) {
				// continue on
			}
		}
		
		if (!f.exists()||!f.canRead()) {
			throw new EncyclopediaException("Can't read file "+f.getAbsolutePath());
		}
		
		// otherwise try readers in order
		for (StripeFileReaderInterface reader : readers) {
			if (reader.canTryToReadFile(f)) {
				return reader.readStripeFile(f, parameters);
			}
		}
		
		throw new EncyclopediaException("Can't read file type "+f.getAbsolutePath());
	}
	
	public static FilenameFilter getFilenameFilter() {
		return new CompoundFilenameFilter(readers);
	}
}
