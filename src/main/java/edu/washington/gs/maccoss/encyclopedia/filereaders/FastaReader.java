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

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class FastaReader {
	public static ArrayList<FastaPeptideEntry> readPeptideFasta(File f) {
		ArrayList<FastaEntryInterface> fasta=readFasta(f);
		
		ArrayList<FastaPeptideEntry> peptides=new ArrayList<FastaPeptideEntry>();
		for (FastaEntryInterface sequence : fasta) {
			peptides.add(sequence.getEntryAsPeptide());
		}
		return peptides;
	}
	
	public static ArrayList<FastaEntryInterface> readFasta(File f) {
		BufferedReader in=null;
		ArrayList<FastaEntryInterface> entryList=new ArrayList<FastaEntryInterface>();
		try {
			in=new BufferedReader(new FileReader(f));
			return readFasta(in, f.getName());

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
	
	public static ArrayList<FastaEntryInterface> readFasta(String s, String fileName) {
		return readFasta(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(s.getBytes(StandardCharsets.UTF_8)))), fileName);
	}
	
	public static ArrayList<FastaEntryInterface> readFasta(InputStream s, String fileName) {
		return readFasta(new BufferedReader(new InputStreamReader(s)), fileName);
	}
	
	public static ArrayList<FastaEntryInterface> readFasta(BufferedReader in, String fileName) {
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
						entryList.add(new FastaEntry(fileName, annotation, sequence.toString()));
					}
					annotation=eachline;
					sequence.setLength(0);
				} else {
					sequence.append(eachline);
				}
			}
			if (annotation!=null) {
				entryList.add(new FastaEntry(fileName, annotation, sequence.toString()));
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
