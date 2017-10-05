package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.UnsupportedEncodingException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.io.OutputMessage;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ExternalExecutor;

public class PercolatorExecutor extends ExternalExecutor {
	public static final String PI_0_TAG="pi_0=";
	public static final String V3_01="v3-01";
	public static final String V2_10="v2-10";
	public static final byte DEFAULT_VERSION_NUMBER=3;

	PercolatorExecutor(File tsv, File outputFile, File decoyFile, int percolatorVersionNumber, boolean useXML) {
		super(generateCommand(tsv, outputFile, decoyFile, percolatorVersionNumber, useXML));
	}

	public static Pair<ArrayList<PercolatorPeptide>, Float> executePercolatorTSV(int percolatorVersionNumber, File featureFile, File percolatorResultFile, File percolatorDecoyFile, float threshold) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		PercolatorExecutor e=new PercolatorExecutor(featureFile, percolatorResultFile, percolatorDecoyFile, percolatorVersionNumber, false);
		BlockingQueue<OutputMessage> result=e.start();

		Float pi0=null;
		String errorMessage=null;
		while (!e.isFinished()||!result.isEmpty()) {
			if (!result.isEmpty()) {
				OutputMessage data=result.take();
				if (!data.isStdOutput()) {
					Logger.logLine(data.getMessage());
					String trim=data.getMessage().trim();
					if (trim.startsWith("Error : ")) {
						errorMessage=trim.substring(8);
					} else if (trim.startsWith("Exception caught: ")) {
						errorMessage=trim.substring(26);
					} else if (trim.indexOf("bad allocation")>=0) {
						errorMessage=trim;
					} else if (trim.startsWith("Selecting pi_0=")) {
						try {
							pi0=Float.parseFloat(trim.substring("Selecting pi_0=".length()));
						} catch (NumberFormatException nfe) {
							Logger.errorLine("Error parsing pi0 from ["+trim+"]");
						}
					}
				}
			} else {
				Thread.sleep(10);
			}
		}

		if (errorMessage!=null) {
			throw new EncyclopediaException(errorMessage);
		}

		checkResult(e);
		
		try {
		    Files.write(percolatorResultFile.toPath(), (PI_0_TAG+pi0+System.lineSeparator()).getBytes(), StandardOpenOption.APPEND);
		}catch (IOException ioe) {
			throw new EncyclopediaException("Error appending to Percolator text file", ioe);
		}

		Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides=PercolatorReader.getPassingPeptidesFromTSV(percolatorResultFile, threshold, false);
		
		return passingPeptides;
	}

	private static void checkResult(PercolatorExecutor e) throws EncyclopediaException {
		if (0 != e.getResultCode()) {
			throw new EncyclopediaException("Percolator exited with non-zero status: " + e.getResultCode());
		}
	}

	static String parsePeptideSequence(String peptideString) {
		return peptideString.substring(peptideString.indexOf('.')+1, peptideString.lastIndexOf('.'));
	}

	static String[] generateCommand(File tsv, File outputFile, File decoyFile, int percolatorVersionNumber, boolean useXML) {
		File percolator=getPercolator(percolatorVersionNumber);

		if (percolatorVersionNumber==2) {
			if (useXML) {
				return new String[] {percolator.getAbsolutePath(), "-y", "--xmloutput", outputFile.getAbsolutePath(), "--decoy-xml-output", tsv.getAbsolutePath()};
			} else {
				return new String[] {percolator.getAbsolutePath(), "--results-peptides", outputFile.getAbsolutePath(), "--decoy-results-peptides", decoyFile.getAbsolutePath(), "-y", tsv.getAbsolutePath()};
			}
		} else {
			if (useXML) {
				return new String[] {percolator.getAbsolutePath(), "-y", "--no-terminate", "-N", "200000", "--xmloutput", outputFile.getAbsolutePath(), "--decoy-xml-output", tsv.getAbsolutePath()};
			} else {
				return new String[] {percolator.getAbsolutePath(), "--results-peptides", outputFile.getAbsolutePath(), "--decoy-results-peptides", decoyFile.getAbsolutePath(), "-y", "--no-terminate", "-N", "200000", tsv.getAbsolutePath()};
			}
		}
	}

	static File getPercolator(int percolatorVersionNumber) {
		String percolatorVersion;
		if (percolatorVersionNumber==2) {
			percolatorVersion=V2_10;
		} else {
			percolatorVersion=V3_01;
		}
		
		try {
			File percolator=File.createTempFile("Percolator-" + percolatorVersion + "-", ".exe");
			percolator.deleteOnExit();
			
			OS os=OSDetector.getOS();
			switch (os) {
				case WINDOWS: {
					InputStream is=PercolatorExecutor.class.getResourceAsStream("/bin/percolator-"+percolatorVersion+".exe");
					Files.copy(is, percolator.toPath(), StandardCopyOption.REPLACE_EXISTING);
					percolator.setExecutable(true);
					
					// not necessary for the crux version of percolator
					//loadLibraryFile(percolator, "xerces-c_3_1.dll");
					//loadLibraryFile(percolator, "msvcr120.dll");
					//loadLibraryFile(percolator, "msvcp120.dll");
					
					return percolator;
				}
				case MAC: {
					InputStream is=PercolatorExecutor.class.getResourceAsStream("/bin/percolator-"+percolatorVersion+".mac");
					Files.copy(is, percolator.toPath(), StandardCopyOption.REPLACE_EXISTING);
					percolator.setExecutable(true);
					return percolator;
				}
				case LINUX:
					InputStream is=PercolatorExecutor.class.getResourceAsStream("/bin/percolator-"+percolatorVersion+".lin");
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