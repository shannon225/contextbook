package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.ExtendedFastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class FastaReader {
	
	public static ArrayList<FastaPeptideEntry> readPeptideFasta(File f, SearchParameters parameters) {
		ArrayList<FastaEntryInterface> fasta=readFasta(f, parameters);
		
		ArrayList<FastaPeptideEntry> peptides=new ArrayList<FastaPeptideEntry>();
		for (FastaEntryInterface sequence : fasta) {
			peptides.add(sequence.getEntryAsPeptide());
		}
		return peptides;
	}

	public static ArrayList<FastaEntryInterface> readFasta(File f, SearchParameters parameters) {
		return readFasta(f, null, parameters);
	}
	public static ArrayList<FastaEntryInterface> readFasta(File f, String keyword, SearchParameters parameters) {
		BufferedReader in=null;
		
		ArrayList<FastaEntryInterface> entryList=new ArrayList<FastaEntryInterface>();
		try {
			in=new BufferedReader(new FileReader(f));
			return readFasta(in, f.getName(), keyword, parameters);

		} catch (IOException ioe) {
			Logger.errorLine("I/O Error found reading FASTA ["+f.getAbsolutePath()+"]");
			Logger.errorException(ioe);
			return entryList;
		} finally {
			if (in!=null) {
				try {
					in.close();
				} catch (IOException ioe) {
					ioe.printStackTrace();
				}
			}
		}
	}
	
	public static ArrayList<FastaEntryInterface> readFasta(String s, String fileName, SearchParameters parameters) {
		return readFasta(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(s.getBytes(StandardCharsets.UTF_8)))), fileName, null, parameters);
	}
	
	public static ArrayList<FastaEntryInterface> readFasta(InputStream s, String fileName, SearchParameters parameters) {
		return readFasta(new BufferedReader(new InputStreamReader(s)), fileName, null, parameters);
	}
	
	//@MoMo this should be removed after adding peffFormat into PARAMETER
	public static ArrayList<FastaEntryInterface> readFasta(BufferedReader in, String fileName, String keyword, SearchParameters parameters) {
		return readFasta(in, fileName, keyword, fileName.toLowerCase().endsWith(".peff"), parameters);
	}
	
	public static ArrayList<FastaEntryInterface> readFasta(BufferedReader in, String fileName, String keyword, boolean peffFormat, SearchParameters parameters) {
		if (peffFormat) {
			Logger.logLine("Reading PEFF formated FASTA database");
		} else {
			Logger.logLine("Reading standard formated FASTA database");
		}
		if (keyword!=null) {
			keyword=keyword.toLowerCase();
		}
		ArrayList<FastaEntryInterface> entryList=new ArrayList<FastaEntryInterface>();
		try {
			String eachline;
			String annotation=null;
			StringBuilder sequence=new StringBuilder();
			while ((eachline=in.readLine())!=null) {
				if (eachline.trim().length()==0) {
					continue;
				}
				if (eachline.startsWith(">")) {
					if (annotation!=null) {
						if (keyword==null||annotation.toLowerCase().indexOf(keyword)>=0) {
							if (peffFormat) {
								entryList.add(new ExtendedFastaEntry(fileName, annotation, sequence.toString(), parameters));	
							}else{
								entryList.add(new FastaEntry(fileName, annotation, sequence.toString()));
							}
						}
					}
					annotation=eachline;
					sequence.setLength(0);
				} else {
					sequence.append(eachline);
				}
			}
			if (annotation!=null) {
				if (keyword==null||annotation.toLowerCase().indexOf(keyword)>=0) {
					if (peffFormat){
						entryList.add(new ExtendedFastaEntry(fileName, annotation, sequence.toString(), parameters));	
					}else{
						entryList.add(new FastaEntry(fileName, annotation, sequence.toString()));
					}
				}
			}
			return entryList;

		} catch (IOException ioe) {
			Logger.errorLine("I/O Error found reading FASTA ["+fileName+"]");
			Logger.errorException(ioe);
			return entryList;
		} finally {
			if (in!=null) {
				try {
					in.close();
				} catch (IOException ioe) {
					ioe.printStackTrace();
				}
			}
		}
	}
}
