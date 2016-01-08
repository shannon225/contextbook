package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.StringTokenizer;
import java.util.concurrent.BlockingQueue;

import com.google.common.base.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;
import edu.washington.gs.maccoss.encyclopedia.utils.io.OutputMessage;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ExternalExecutor;

public class PercolatorExecutor extends ExternalExecutor {

	PercolatorExecutor(File tsv) {
		super(generateCommand(tsv, Optional.fromNullable((File)null)));
	}
	PercolatorExecutor(File tsv, Optional<File> percolatorLocation) {
		super(generateCommand(tsv, percolatorLocation));
	}
	
	public static ArrayList<ScoredObject<String>> executePercolator(Optional<File> percolatorLocation, File featureFile, File outputFile, float threshold) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		PercolatorExecutor e=new PercolatorExecutor(featureFile, percolatorLocation);
		BlockingQueue<OutputMessage> result=e.start();
		
		boolean isFirst=true;
		boolean record=true;
		PrintWriter writer=new PrintWriter(outputFile, "UTF-8");
		ArrayList<ScoredObject<String>> passingPeptides=new ArrayList<ScoredObject<String>>();
		while (!e.isFinished()||!result.isEmpty()) {
			if (!result.isEmpty()) {
				OutputMessage data=result.take();
				if (data.isStdOutput()) {
					if (isFirst) {
						isFirst=false;
					} else if (record) {
						StringTokenizer st=new StringTokenizer(data.getMessage());
						String psmID=st.nextToken(); // PSMid
						st.nextToken(); // score
						float qvalue=Float.parseFloat(st.nextToken()); //Q-value
						//st.nextToken(); // PEP
						//String peptideString=st.nextToken();
						//String peptideSequence = parsePeptideSequence(peptideString);
						
						if (qvalue<threshold) {
							ScoredObject<String> peptide=new ScoredObject<String>(qvalue, psmID);
							passingPeptides.add(peptide);
						} else {
							record=false;
						}
					}
					writer.println(data.getMessage());
				} else {
					Logger.logLine(data.getMessage());
				}
			} else {
				Thread.sleep(10);
			}
		}
		writer.flush();
		writer.close();
		return passingPeptides;
	}

	static String parsePeptideSequence(String peptideString) {
		return peptideString.substring(peptideString.indexOf('.')+1, peptideString.lastIndexOf('.'));
	}
	
	static String[] generateCommand(File tsv, Optional<File> percolatorLocation) {
		File percolator=getPercolator(percolatorLocation);
		
		return new String[] {percolator.getAbsolutePath(), tsv.getAbsolutePath()};
	}
	
	static File getPercolator(Optional<File> percolatorLocation) {
		if (percolatorLocation.isPresent()) return percolatorLocation.get();
		
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
					InputStream is=PercolatorExecutor.class.getResourceAsStream("/bin/percolator-v2-08.lin");
					Files.copy(is, percolator.toPath(), StandardCopyOption.REPLACE_EXISTING);
					percolator.setExecutable(true);
					return percolator;
			}
			throw new EncyclopediaException("Sorry, Percolator for "+OSDetector.getOSName(os)+" is not set up yet!");
		} catch (IOException ioe) {
			throw new EncyclopediaException("Unexpected exception finding Percolator", ioe);
		}
	}

	static void loadLibraryFile(File percolator, String target) throws IOException {
		File file=new File(percolator.getParentFile(), target);
		file.deleteOnExit();
		InputStream is=PercolatorExecutor.class.getResourceAsStream("/bin/"+target);
		Files.copy(is, file.toPath(), StandardCopyOption.REPLACE_EXISTING);
	}
}