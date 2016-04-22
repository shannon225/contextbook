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
import java.util.Optional;
import java.util.StringTokenizer;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;
import edu.washington.gs.maccoss.encyclopedia.utils.io.OutputMessage;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ExternalExecutor;

public class PercolatorExecutor extends ExternalExecutor {

	PercolatorExecutor(File tsv, File outputFile, boolean useXML) {
		super(generateCommand(tsv, outputFile, Optional.ofNullable((File)null), useXML));
	}
	PercolatorExecutor(File tsv, File outputFile, Optional<File> percolatorLocation, boolean useXML) {
		super(generateCommand(tsv, outputFile, percolatorLocation, useXML));
	}
	
	public static ArrayList<ScoredObject<String>> executePercolatorXML(Optional<File> percolatorLocation, File featureFile, File percolatorResultFile, float threshold) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		PercolatorExecutor e=new PercolatorExecutor(featureFile, percolatorResultFile, percolatorLocation, true);
		BlockingQueue<OutputMessage> result=e.start();
		
		while (!e.isFinished()||!result.isEmpty()) {
			if (!result.isEmpty()) {
				OutputMessage data=result.take();
				if (!data.isStdOutput()) {
					Logger.logLine(data.getMessage());
				}
			} else {
				Thread.sleep(10);
			}
		}

		ArrayList<ScoredObject<String>> passingPeptides=PercolatorReader.getPassingPeptidesFromXML(percolatorResultFile, threshold);
		
		return passingPeptides;
	}
	
	public static ArrayList<ScoredObject<String>> executePercolatorTSV(Optional<File> percolatorLocation, File featureFile, File outputFile, float threshold) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		PercolatorExecutor e=new PercolatorExecutor(featureFile, outputFile, percolatorLocation, false);
		BlockingQueue<OutputMessage> result=e.start();
		
		String errorMessage=null;
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
					String trim=data.getMessage().trim();
					if (trim.startsWith("Error : ")) {
						errorMessage=trim.substring(8);
					}
					if (trim.startsWith("bad allocation")) {
						errorMessage=trim;
					}
				}
			} else {
				Thread.sleep(10);
			}
		}
		writer.flush();
		writer.close();
		if (errorMessage!=null) {
			throw new EncyclopediaException(errorMessage);
		}
		return passingPeptides;
	}

	static String parsePeptideSequence(String peptideString) {
		return peptideString.substring(peptideString.indexOf('.')+1, peptideString.lastIndexOf('.'));
	}
	
	static String[] generateCommand(File tsv, File outputFile, Optional<File> percolatorLocation, boolean useXML) {
		File percolator=getPercolator(percolatorLocation);

		OS os=OSDetector.getOS();
		switch (os) {
			case MAC: {
				if (useXML) {
					return new String[] {percolator.getAbsolutePath(), "-X", outputFile.getAbsolutePath(), "--decoy-xml-output", tsv.getAbsolutePath()};
				} else {
					return new String[] {percolator.getAbsolutePath(), tsv.getAbsolutePath()};
				}
			}
			default: {
				if (useXML) {
					return new String[] {percolator.getAbsolutePath(), "-y", "-X", outputFile.getAbsolutePath(), "--decoy-xml-output", tsv.getAbsolutePath()};
				} else {
					return new String[] {percolator.getAbsolutePath(), "-y", tsv.getAbsolutePath()};
				}	
			}
		}
	}
	
	static File getPercolator(Optional<File> percolatorLocation) {
		if (percolatorLocation.isPresent()) return percolatorLocation.get();
		
		try {
			File percolator=File.createTempFile("Percolator", ".exe");
			percolator.deleteOnExit();
			
			OS os=OSDetector.getOS();
			switch (os) {
				case WINDOWS: {
					InputStream is=PercolatorExecutor.class.getResourceAsStream("/bin/percolator-v2-10.exe");
					Files.copy(is, percolator.toPath(), StandardCopyOption.REPLACE_EXISTING);
					percolator.setExecutable(true);
					
					// not necessary for the crux version of percolator
					//loadLibraryFile(percolator, "xerces-c_3_1.dll");
					//loadLibraryFile(percolator, "msvcr120.dll");
					//loadLibraryFile(percolator, "msvcp120.dll");
					
					return percolator;
				}
				case MAC: {
					InputStream is=PercolatorExecutor.class.getResourceAsStream("/bin/percolator-v2-06.mac");
					Files.copy(is, percolator.toPath(), StandardCopyOption.REPLACE_EXISTING);
					percolator.setExecutable(true);
					return percolator;
				}
				case LINUX:
					InputStream is=PercolatorExecutor.class.getResourceAsStream("/bin/percolator-v2-10.lin");
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