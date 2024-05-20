package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.FilenameFilter;
import java.lang.ref.WeakReference;
import java.util.HashMap;
import java.util.Optional;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.gui.general.CompoundFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class StripeFileGenerator {

	private static final HashMap<File, WeakReference<StripeFileInterface>> loadedFiles=new HashMap<>();
	
	private static final DIAFileReader DIA_FILE_READER=new DIAFileReader();
	private static final MzmlToDIAConverter MZML_TO_DIA_CONVERTER = new MzmlToDIAConverter();
	private static final StripeFileReaderInterface[] readers=new StripeFileReaderInterface[] {
			DIA_FILE_READER,
			MZML_TO_DIA_CONVERTER,
			new MGFToDIAConverter()
	};

	public static StripeFileInterface getFile(File f, SearchParameters parameters) {
		return getFile(f, parameters, false);
	}

	public static StripeFileInterface getFile(File f, SearchParameters parameters, boolean isOpenFileInPlace) {
		if (loadedFiles.containsKey(f)) {
			StripeFileInterface maybeLoaded=loadedFiles.get(f).get();
			if (maybeLoaded==null||!maybeLoaded.isOpen()) {
				loadedFiles.remove(f);
			} else {
				return maybeLoaded;
			}
		}

		StripeFileInterface file=innerGetFile(f, parameters, isOpenFileInPlace);
		loadedFiles.put(f, new WeakReference<StripeFileInterface>(file));
		return file;
	}
	
	private static StripeFileInterface innerGetFile(File f, SearchParameters parameters, boolean isOpenFileInPlace) {	
		// try to change name to .DIA and read
		String absolutePath=f.getAbsolutePath();
		File diaFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+StripeFile.DIA_EXTENSION);
		if (diaFile.exists()&&diaFile.canRead()) {
			try {
				return DIA_FILE_READER.readStripeFile(diaFile, parameters, isOpenFileInPlace);
			} catch (EncyclopediaException ee) {
				// continue on
			}
		}
		
		// try to change name to .mzML and read
		File foundMZMLFile=null;
		File mzMLFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+".mzml");
		if (mzMLFile.exists()&&mzMLFile.canRead()) {
			foundMZMLFile=mzMLFile;
		}
		if (foundMZMLFile==null) {
			mzMLFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+".mzML");
			if (mzMLFile.exists()&&mzMLFile.canRead()) {
				foundMZMLFile=mzMLFile;
			}
		}
		if (foundMZMLFile==null) {
			mzMLFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+".MZML");
			if (mzMLFile.exists()&&mzMLFile.canRead()) {
				foundMZMLFile=mzMLFile;
			}
		}
		if (foundMZMLFile!=null&&foundMZMLFile.exists()&&foundMZMLFile.canRead()) {
			try {
				return MZML_TO_DIA_CONVERTER.readStripeFile(foundMZMLFile, parameters, isOpenFileInPlace);
			} catch (EncyclopediaException ee) {
				// continue on
			}
		}
		
		Optional<String> potentialName=getBuggyFileName(absolutePath);
		if (potentialName.isPresent()) {
			try {
				String alternateName=potentialName.get();
				Logger.logLine("Found alternative raw file name, considering: "+alternateName);
				return getFile(new File(alternateName), parameters);
			} catch (EncyclopediaException ee) {
				// ignore in case this isn't actually a buggy file name
			}
		}
		
		if (!f.exists()||!f.canRead()) {
			throw new EncyclopediaException("Can't read file "+f.getAbsolutePath());
		}
		
		// otherwise try readers in order
		for (StripeFileReaderInterface reader : readers) {
			if (reader.canTryToReadFile(f)) {
				return reader.readStripeFile(f, parameters, isOpenFileInPlace);
			}
		}
		
		throw new EncyclopediaException("Can't read file type "+f.getAbsolutePath());
	}
	
	static Optional<String> getBuggyFileName(String f) {
		Pattern p=Pattern.compile("^(.*\\.[mM][zZ][mM][lL])[0-9]*\\.mzml$");

		Matcher m=p.matcher(f);
		if (m.find()) {
			return Optional.of(m.group(1));
		} else {
			return Optional.empty();
		}
	}
	
	public static FilenameFilter getFilenameFilter() {
		return new CompoundFilenameFilter(readers);
	}
}
