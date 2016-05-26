package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.Version;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class BlibToLibraryConverter {

	public static LibraryInterface getFile(File f) {
		if (!f.exists()||!f.canRead()) {
			throw new EncyclopediaException("Can't read file "+f.getAbsolutePath());
		}
		
		// first try to read if .ELIB
		if (f.getName().toLowerCase().endsWith(LibraryFile.ELIB)) {
			Optional<LibraryInterface> optional=openLibraryFile(f);
			if (optional.isPresent()) return optional.get();
		}
		
		// then try to change name to .ELIB and read
		String absolutePath=f.getAbsolutePath();
		File libraryFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.ELIB);
		if (libraryFile.exists()&&libraryFile.canRead()) {
			Optional<LibraryInterface> optional=openLibraryFile(libraryFile);
			if (optional.isPresent()) return optional.get();
		}

		throw new EncyclopediaException("Can't read file type "+f.getAbsolutePath());
	}

	static Optional<LibraryInterface> openLibraryFile(File f) {
		try {
			LibraryFile elibFile=new LibraryFile();
			elibFile.openFile(f);
			Version version=elibFile.getVersion();
			if (LibraryFile.isVersionAcceptable(version)) {
				return Optional.of((LibraryInterface)elibFile);
			} else {
				Logger.errorLine(f.getName()+" is an old Library file (version:"+version+")! Please delete it!");
			}
			return Optional.empty();
			
		} catch (IOException ioe) {
			throw new EncyclopediaException("Error reading ELIB file!", ioe);
		} catch (SQLException sqle) {
			throw new EncyclopediaException("Error reading ELIB file!", sqle);
		}
	}

	public static LibraryInterface convert(File blibFile, Optional<File> irtFile, File fastaFile) {
		String absolutePath=blibFile.getAbsolutePath();
		File elibFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.ELIB);
		return convert(blibFile, elibFile, irtFile, fastaFile);
	}
	
	static LibraryInterface convert(File blibFile, File elibFile, Optional<File> irtFile, File fastaFile) {
		TObjectFloatHashMap<String> irtMap=null;
		try {
			Logger.logLine("Indexing "+blibFile.getName()+" ...");
			LibraryFile elib=new LibraryFile();
			elib.openFile();
			
			BlibFile blib = new BlibFile();
			blib.openFile(blibFile);
			if (irtFile.isPresent()) {
				IRTdbFile irt=new IRTdbFile(irtFile.get());
				irtMap=irt.getIRTs();
			}
			blib.getCopyEntriesToLibrary(elib, Optional.ofNullable(irtMap), fastaFile);	
			elib.saveAsFile(elibFile);
			return elib;
			
		} catch (DataFormatException dfe) {
			throw new EncyclopediaException("ELIB writing Data Formatting error!", dfe);
		} catch (IOException ioe) {
			throw new EncyclopediaException("ELIB writing IO error!", ioe);
		} catch (SQLException sqle) {
			throw new EncyclopediaException("ELIB writing SQL error!", sqle);
		}
	}
}
