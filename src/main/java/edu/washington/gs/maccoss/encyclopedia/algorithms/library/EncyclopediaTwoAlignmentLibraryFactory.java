package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.HashMap;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;

public class EncyclopediaTwoAlignmentLibraryFactory {
	// echos DigestionEnzyme
	private static final HashMap<String, String> enzymeNames=new HashMap<String, String>();
	private static void prepopulateNames() {
		enzymeNames.put("Trypsin", "trimmed_Trypsin.csv.dlib");
		enzymeNames.put("Trypsin/p", "trimmed_Trypsin.csv.dlib");
		enzymeNames.put("Lys-C", "trimmed_Lys-C.csv.dlib");
		enzymeNames.put("Lys-N", "trimmed_Lys-N.csv.dlib");
		enzymeNames.put("Arg-C", "trimmed_Arg-C.csv.dlib");
		enzymeNames.put("Glu-C", "trimmed_Glu-C.csv.dlib");
		enzymeNames.put("Chymotrypsin", "trimmed_Chymotrypsin.csv.dlib");
		enzymeNames.put("Pepsin A", "trimmed_Pepsin.csv.dlib");
		enzymeNames.put("Elastase", "trimmed_Elastase.csv.dlib");
		enzymeNames.put("Thermolysin", "trimmed_Thermolysin.csv.dlib");
		
	}
	
	public static Optional<File> getPreAlignmentFile(DigestionEnzyme enzyme) {
		if (enzymeNames.size()==0) prepopulateNames();
		
		String fname=enzymeNames.get(enzyme.getName());
		Logger.logLine("Looking for pre-alignment library for ["+enzyme.getName()+"]. Found ["+fname+"]");
		
		if (fname==null) {
			return Optional.empty();
		} else {
			try {
				File library = File.createTempFile(fname, ".dlib");
				library.deleteOnExit();
				InputStream is = PercolatorExecutor.class.getResourceAsStream("/libraries/"+fname);
				Files.copy(is, library.toPath(), StandardCopyOption.REPLACE_EXISTING);
				return Optional.of(library);
				
			} catch (IOException ioe) {
				throw new EncyclopediaException("Unexpected exception finding pre-alignment library", ioe);
			} 
		}
	}
}
