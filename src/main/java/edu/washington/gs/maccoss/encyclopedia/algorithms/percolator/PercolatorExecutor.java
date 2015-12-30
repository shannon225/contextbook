package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ExternalExecutor;

public class PercolatorExecutor extends ExternalExecutor {
	
	public PercolatorExecutor(File tsv) {
		super(generateCommand(tsv));
	}
	
	public static String[] generateCommand(File tsv) {
		File percolator=getPercolator();
		
		return new String[] {percolator.getAbsolutePath(), tsv.getAbsolutePath()};
	}
	
	public static File getPercolator() {
		try {
			File percolator=File.createTempFile("Percolator", ".exe");
			percolator.deleteOnExit();
			
			OS os=OSDetector.getOS();
			switch (os) {
				case WINDOWS: {
					InputStream is=PercolatorExecutor.class.getResourceAsStream("/bin/percolator-v2-06.exe");
					Files.copy(is, percolator.toPath(), StandardCopyOption.REPLACE_EXISTING);
					percolator.setExecutable(true);
					
					loadLibraryFile(percolator, "xerces-c_3_1.dll");
					loadLibraryFile(percolator, "msvcr120.dll");
					loadLibraryFile(percolator, "msvcp120.dll");
					
					return percolator;
				}
				case MAC: {
					InputStream is=PercolatorExecutor.class.getResourceAsStream("/bin/percolator-v2-06.mac");
					Files.copy(is, percolator.toPath(), StandardCopyOption.REPLACE_EXISTING);
					percolator.setExecutable(true);
					return percolator;
				}
				case LINUX:
					InputStream is=PercolatorExecutor.class.getResourceAsStream("/bin/percolator-v2-04.lin");
					Files.copy(is, percolator.toPath(), StandardCopyOption.REPLACE_EXISTING);
					percolator.setExecutable(true);
					return percolator;
			}
			throw new EncyclopediaException("Sorry, Percolator for "+OSDetector.getOSName(os)+" is not set up yet!");
		} catch (IOException ioe) {
			throw new EncyclopediaException("Unexpected exception finding Percolator", ioe);
		}
	}

	public static void loadLibraryFile(File percolator, String target) throws IOException {
		File file=new File(percolator.getParentFile(), target);
		file.deleteOnExit();
		InputStream is=PercolatorExecutor.class.getResourceAsStream("/bin/"+target);
		Files.copy(is, file.toPath(), StandardCopyOption.REPLACE_EXISTING);
	}
}