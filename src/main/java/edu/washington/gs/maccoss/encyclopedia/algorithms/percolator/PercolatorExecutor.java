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

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;
import edu.washington.gs.maccoss.encyclopedia.utils.io.OutputMessage;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ExternalExecutor;

public class PercolatorExecutor extends ExternalExecutor {
	public static final String V3_01="v3-01";
	public static final String V2_10="v2-10";
	public static final byte DEFAULT_VERSION_NUMBER=3;

	PercolatorExecutor(File tsv, File outputFile, int percolatorVersionNumber, boolean useXML) {
		super(generateCommand(tsv, outputFile, percolatorVersionNumber, useXML));
	}

	public static ArrayList<ScoredObject<String>> executePercolatorXML(int percolatorVersionNumber, File featureFile, File percolatorResultFile, float threshold) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		PercolatorExecutor e=new PercolatorExecutor(featureFile, percolatorResultFile, percolatorVersionNumber, true);
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

		checkResult(e);

		ArrayList<ScoredObject<String>> passingPeptides=PercolatorReader.getPassingPeptidesFromXML(percolatorResultFile, threshold);
		
		return passingPeptides;
	}
	
	public static ArrayList<PercolatorPeptide> executePercolatorTSV(int percolatorVersionNumber, File featureFile, File outputFile, float threshold) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		PercolatorExecutor e=new PercolatorExecutor(featureFile, outputFile, percolatorVersionNumber, false);
		BlockingQueue<OutputMessage> result=e.start();
		
		String errorMessage=null;
		boolean isFirst=true;
		boolean record=true;
		PrintWriter writer=new PrintWriter(outputFile, "UTF-8");
		ArrayList<PercolatorPeptide> passingPeptides=new ArrayList<PercolatorPeptide>();
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
						float pep=Float.parseFloat(st.nextToken()); // PEP
						st.nextToken(); // peptide
						String proteinIds=st.nextToken();
						//String peptideString=st.nextToken();
						//String peptideSequence = parsePeptideSequence(peptideString);
						
						if (qvalue<threshold) {
							if (!proteinIds.startsWith(LibraryEntry.DECOY_STRING)) {
								passingPeptides.add(new PercolatorPeptide(psmID, proteinIds, qvalue, pep));
							}
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
					} else if (trim.startsWith("Exception caught: ")) {
						errorMessage=trim.substring(26);
					} else if (trim.indexOf("bad allocation")>=0) {
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

		checkResult(e);

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

	@SuppressWarnings("unused")
	static String[] generateCommand(File tsv, File outputFile, int percolatorVersionNumber, boolean useXML) {
		File percolator=getPercolator(percolatorVersionNumber);

		if (percolatorVersionNumber==2) {
			if (useXML) {
				return new String[] {percolator.getAbsolutePath(), "-y", "-X", outputFile.getAbsolutePath(), "--decoy-xml-output", tsv.getAbsolutePath()};
			} else {
				return new String[] {percolator.getAbsolutePath(), "-y", tsv.getAbsolutePath()};
			}
		} else {
			if (useXML) {
				return new String[] {percolator.getAbsolutePath(), "-y", "--no-terminate", "-N", "200000", "-X", outputFile.getAbsolutePath(), "--decoy-xml-output", tsv.getAbsolutePath()};
			} else {
				return new String[] {percolator.getAbsolutePath(), "-y", "--no-terminate", "-N", "200000", tsv.getAbsolutePath()};
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