package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.Version;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class BlibToLibraryConverter {

	/**
	 * Utility for opening an existing ELIB/DLIB file, located by checking
	 * if the given file is of an appropriate format or if a file with the
	 * same name but the appropriate extension exists. Does not convert BLIBs!
	 * <p>
	 * TODO: consider moving this utility to a more appropriately-named class
	 *
	 * @param f The file for which an existing ELIB/DLIB should be located.
	 * @return An {@code LibraryInterface} corresponding to the given file.
	 * Note that the returned instance may be backed by another file
	 * than the given one.
	 * @throws EncyclopediaException if an appropriate file can't be located,
	 *                               or if one can but can't be used.
	 */
	public static LibraryInterface getFile(File f) {
		if (!f.exists() || !f.canRead()) {
			throw new EncyclopediaException("Can't read file " + f.getAbsolutePath());
		}

		// first try to read if .ELIB
		if (f.getName().toLowerCase().endsWith(LibraryFile.ELIB)) {
			Optional<? extends LibraryInterface> optional = openLibraryFile(f);
			if (optional.isPresent()) return optional.get();
		}

		// then try to change name to .ELIB and read
		String absolutePath = f.getAbsolutePath();
		File libraryFile = new File(absolutePath.substring(0, absolutePath.lastIndexOf('.')) + LibraryFile.ELIB);
		if (libraryFile.exists() && libraryFile.canRead()) {
			Optional<? extends LibraryInterface> optional = openLibraryFile(libraryFile);
			if (optional.isPresent()) return optional.get();
		}

		// try to read if .DLIB
		if (f.getName().toLowerCase().endsWith(LibraryFile.DLIB)) {
			Optional<? extends LibraryInterface> optional = openLibraryFile(f);
			if (optional.isPresent()) return optional.get();
		}

		// then try to change name to .DLIB and read
		libraryFile = new File(absolutePath.substring(0, absolutePath.lastIndexOf('.')) + LibraryFile.DLIB);
		if (libraryFile.exists() && libraryFile.canRead()) {
			Optional<? extends LibraryInterface> optional = openLibraryFile(libraryFile);
			if (optional.isPresent()) return optional.get();
		}

		throw new EncyclopediaException("Can't read file type " + f.getAbsolutePath());
	}

	static Optional<LibraryFile> openLibraryFile(File f) {
		try {
			LibraryFile elibFile = new LibraryFile();
			elibFile.openFile(f);
			Version version = elibFile.getVersion();
			if (LibraryFile.isVersionAcceptable(version)) {
				return Optional.of(elibFile);
			} else {
				if (version.amIAbove(LibraryFile.MOST_RECENT_VERSION)) {
					Logger.errorLine("The library file " + f.getName() + " is newer than expected (version:" + version + ")! Please download a newer upgrade that supports this library!");
				} else {
					Logger.errorLine("The library file " + f.getName() + " is too old to be used (version:" + version + ")! Please delete it!");
				}
			}
			return Optional.empty();

		} catch (IOException ioe) {
			throw new EncyclopediaException("Error reading ELIB file!", ioe);
		} catch (SQLException sqle) {
			throw new EncyclopediaException("Error reading ELIB file!", sqle);
		}
	}

	public static LibraryInterface convert(File blibFile, Optional<File> irtFile, File fastaFile, SearchParameters params) {
		String absolutePath = blibFile.getAbsolutePath();
		File elibFile = new File(absolutePath.substring(0, absolutePath.lastIndexOf('.')) + LibraryFile.DLIB);
		return convert(blibFile, elibFile, irtFile, fastaFile, params);
	}

	static LibraryInterface convert(File blibFile, File elibFile, Optional<File> irtFile, File fastaFile, SearchParameters params) {
		TObjectFloatHashMap<String> irtMap = null;
		try {
			Logger.logLine("Indexing " + blibFile.getName() + " ...");
			LibraryFile elib = new LibraryFile();
			elib.openFile();

			BlibFile blib = new BlibFile();
			blib.openFile(blibFile);
			if (irtFile.isPresent()) {
				IRTdbFile irt = new IRTdbFile(irtFile.get());
				irtMap = irt.getIRTs();
			}
			blib.getCopyEntriesToLibrary(elib, Optional.ofNullable(irtMap), fastaFile, params);
			elib.saveAsFile(elibFile);
			elib.close();
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
