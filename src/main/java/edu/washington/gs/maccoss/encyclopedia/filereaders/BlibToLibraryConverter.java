package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.Version;
import gnu.trove.map.hash.TObjectFloatHashMap;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Optional;
import java.util.zip.DataFormatException;

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
		return internalGetLibraryFile(f);
	}
	
	public static LibraryInterface getFile(File f, File fastaFile, SearchParameters parameters) {
		LibraryFile libraryFile=internalGetLibraryFile(f);
		try {
			if (LibraryEntryCleaner.doesLibraryNeedCleaning(libraryFile)) {
				String name = f.getName();
				File newLibrary=new File(f.getParentFile(), "cleaned_"+name.substring(0, name.lastIndexOf('.'))+LibraryFile.DLIB);
				if (newLibrary.exists()) {
					return getFile(newLibrary);
				} else {
					libraryFile=LibraryEntryCleaner.cleanLibrary(true, false, f, newLibrary, fastaFile, parameters);
				}
			};
		} catch (SQLException sqle) {
		} catch (IOException ioe) {
		} catch (DataFormatException dfe) {
		}
		
		return libraryFile;
	}

	private static LibraryFile internalGetLibraryFile(File f) {
		if (!f.exists() || !f.canRead()) {
			throw new EncyclopediaException("Can't read library file " + f.getAbsolutePath());
		}

		// first try to read if .ELIB
		if (f.getName().toLowerCase().endsWith(LibraryFile.ELIB)) {
			Optional<LibraryFile> optional = openLibraryFile(f);
			if (optional.isPresent()) return optional.get();
		}

		// then try to change name to .ELIB and read
		String absolutePath = f.getAbsolutePath();
		File libraryFile = new File(absolutePath.substring(0, absolutePath.lastIndexOf('.')) + LibraryFile.ELIB);
		if (libraryFile.exists() && libraryFile.canRead()) {
			Optional<LibraryFile> optional = openLibraryFile(libraryFile);
			if (optional.isPresent()) return optional.get();
		}

		// try to read if .DLIB
		if (f.getName().toLowerCase().endsWith(LibraryFile.DLIB)) {
			Optional<LibraryFile> optional = openLibraryFile(f);
			if (optional.isPresent()) return optional.get();
		}

		// then try to change name to .DLIB and read
		libraryFile = new File(absolutePath.substring(0, absolutePath.lastIndexOf('.')) + LibraryFile.DLIB);
		if (libraryFile.exists() && libraryFile.canRead()) {
			Optional<LibraryFile> optional = openLibraryFile(libraryFile);
			if (optional.isPresent()) return optional.get();
		}

		throw new EncyclopediaException("Can't read library file type " + f.getAbsolutePath());
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

	/**
	 * @return A <emph>closed</emph> {@code LibraryFile} instance pointing to {@code elibFile}
	 *         containing the results of conversion.
	 */
	public static LibraryFile convert(File blibFile, Optional<File> irtFile, File fastaFile, boolean higherScoreBetter, SearchParameters params) {
		String absolutePath = blibFile.getAbsolutePath();
		File elibFile = new File(absolutePath.substring(0, absolutePath.lastIndexOf('.')) + LibraryFile.DLIB);
		return convert(blibFile, elibFile, irtFile, fastaFile, higherScoreBetter, params);
	}

	/**
	 * @return A <emph>closed</emph> {@code LibraryFile} instance pointing to {@code elibFile}
	 *         containing the results of conversion.
	 */
	public static LibraryFile convert(File blibFile, File elibFile, Optional<File> irtFile, File fastaFile, boolean higherScoreBetter, SearchParameters params) {
		TObjectFloatHashMap<String> irtMap = null;
		try {
			Logger.logLine("Indexing " + blibFile.getName() + " ...");
			LibraryFile elib = new LibraryFile();
			elib.openFile();

			BlibFile blib = new BlibFile();
			blib.openFile(blibFile);
			if (irtFile.isPresent()) {
				IRTdbFile irt = new IRTdbFile();
				irt.openFile(irtFile.get());
				irtMap = irt.getIRTs();
			}
			blib.getCopyEntriesToLibrary(elib, Optional.ofNullable(irtMap), fastaFile, params);

			// if there's no iRT file, then try to adjust the RTs
//			if (!irtFile.isPresent()) {
//				ArrayList<LibraryEntry> entries=elib.getAllEntries(false, params.getAAConstants());
//				entries=LibraryEntryCleaner.correctRTs(entries, elibFile);
//				entries=LibraryEntryCleaner.removeDuplicateEntries(entries, higherScoreBetter);
//
//				LibraryFile elib2 = new LibraryFile();
//				elib2.openFile();
//				elib2.dropIndices();
//				elib2.addEntries(entries);
//				elib2.addProteinsFromEntries(entries);
//				elib2.createIndices();
//				elib=elib2;
//			}
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
