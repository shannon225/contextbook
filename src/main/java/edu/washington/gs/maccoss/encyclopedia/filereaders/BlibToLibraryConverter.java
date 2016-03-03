package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class BlibToLibraryConverter {
	public static void main(String[] args) throws IOException, SQLException, DataFormatException {
		File[] blibFiles;
		File iRTLibraryFile;
		/*
		File blibDir=new File("/Users/searleb/Documents/villen_manuscript/phospho/");
		File[] blibFiles=new File[] {new File(blibDir, "GerberS.blib"), new File(blibDir, "MannM_HumanPhospho.blib"), new File(blibDir, "VillenJ_Velos_1to200.blib"), new File(blibDir, "HeckA.blib"),
				new File(blibDir, "VillenJ_Exactive_HumanPhosphoproteome.blib"), new File(blibDir, "VillenJ_Velos_201to424_HumanPhospho.blib")};

		
		//blibFiles=new File[] {new File("/Users/searleb/Documents/projects/encyclopedia/mzml/121115_BCS_HeLa_24mz_400_1000.blib")};
		//libraryFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/121115_BCS_HeLa_24mz_400_1000.elib");
		 * 
		blibDir=new File("/Users/searleb/Documents/school/projects/");
		blibFiles=new File[] {new File(blibDir, "VillenJ_Exactive_HumanPhosphoproteome.blib")};
		blibFiles=new File[] {new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/VillenJ_Exactive_HumanPhosphoproteome.blib")};
		libraryFile=new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/VillenJ_Exactive_HumanPhosphoproteome.elib");
		iRTLibraryFile=new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/HumanPhosphoProteomeRT.irtdb");
		*/
		
		blibFiles=new File[] {new File("/Users/searleb/Documents/projects/encyclopedia/mzml/momo/20151125_122sample_id_iRT.blib")};
		iRTLibraryFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/momo/20151125_122sample_id_iRT.irtdb");
		getFile(blibFiles[0], Optional.ofNullable(iRTLibraryFile));
	}

	public static LibraryInterface getFile(File f, Optional<File> irtFile) {
		if (!f.exists()||!f.canRead()) {
			throw new EncyclopediaException("Can't read file "+f.getAbsolutePath());
		}
		
		// first try to read if .ELIB
		if (f.getName().toLowerCase().endsWith(LibraryFile.ELIB)) {
			return openLibraryFile(f);
		}
		
		// then try to change name to .ELIB and read
		String absolutePath=f.getAbsolutePath();
		File libraryFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.ELIB);
		if (libraryFile.exists()&&libraryFile.canRead()) {
			return openLibraryFile(libraryFile);
		}
		
		// otherwise check for BLIB and convert
		if (f.getName().toLowerCase().endsWith(BlibFile.BLIB)) {
			return convert(f, libraryFile, irtFile);
		} else {
			throw new EncyclopediaException("Can't read file type "+f.getAbsolutePath());
		}
	}

	static LibraryInterface openLibraryFile(File f) {
		try {
			LibraryFile elibFile=new LibraryFile();
			elibFile.openFile(f);
			return elibFile;
		} catch (IOException ioe) {
			throw new EncyclopediaException("Error reading ELIB file!", ioe);
		} catch (SQLException sqle) {
			throw new EncyclopediaException("Error reading ELIB file!", sqle);
		}
	}

	static LibraryInterface convert(File blibFile, File elibFile, Optional<File> irtFile) {
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
			blib.getStreamEntriesToLibrary(elib, Optional.ofNullable(irtMap));	
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
