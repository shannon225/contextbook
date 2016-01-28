package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.zip.DataFormatException;

import gnu.trove.map.hash.TObjectFloatHashMap;

public class BlibToLibraryConverter {
	public static void main(String[] args) throws IOException, SQLException, DataFormatException {
		File blibDir=new File("/Users/searleb/Documents/villen_manuscript/phospho/");
		File[] blibFiles=new File[] {new File(blibDir, "GerberS.blib"), new File(blibDir, "MannM_HumanPhospho.blib"), new File(blibDir, "VillenJ_Velos_1to200.blib"), new File(blibDir, "HeckA.blib"),
				new File(blibDir, "VillenJ_Exactive_HumanPhosphoproteome.blib"), new File(blibDir, "VillenJ_Velos_201to424_HumanPhospho.blib")};

		blibDir=new File("/Users/searleb/Documents/school/projects/");
		blibFiles=new File[] {new File(blibDir, "VillenJ_Exactive_HumanPhosphoproteome.blib")};
		blibFiles=new File[] {new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/VillenJ_Exactive_HumanPhosphoproteome.blib")};
		//blibFiles=new File[] {new File("/Users/searleb/Documents/projects/encyclopedia/mzml/hela_6mz.blib")};
		
		File libraryFile=new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/VillenJ_Exactive_HumanPhosphoproteome.elib");
		//File libraryFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/hela_6mz.elib");
		
		File iRTLibraryFile=new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/phosphoiRT.irtdb");
		IRTdbFile irt=new IRTdbFile(iRTLibraryFile);
		TObjectFloatHashMap<String> irtMap=irt.getIRTs();

		LibraryFile library=new LibraryFile();
		library.openFile();
		for (File blibFile : blibFiles) {
			BlibFile blib = new BlibFile();
			blib.openFile(blibFile);
			blib.getStreamEntriesToLibrary(library, irtMap);	
		}
		library.saveAsFile(libraryFile);
	}
}
