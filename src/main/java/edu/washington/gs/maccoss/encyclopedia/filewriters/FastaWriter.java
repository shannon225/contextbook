package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;

import org.apache.commons.lang3.text.WordUtils;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class FastaWriter {
	private final PrintWriter writer;
	public FastaWriter(File f) {
		try {
			writer=new PrintWriter(f, "UTF-8");
		} catch (FileNotFoundException|UnsupportedEncodingException e) {
			throw new EncyclopediaException("Error writing FASTA!", e);
		}
	}
	
	public void write(FastaEntryInterface entry) {
		writer.print('>');
		writer.println(entry.getAnnotation());
		writer.println(WordUtils.wrap(entry.getSequence(), 80, System.lineSeparator(), true));
	}
	
	public void close() {
		writer.close();
	}
	
	public static void writeFasta(File f, ArrayList<FastaEntryInterface> list) {
		FastaWriter writer=new FastaWriter(f);
		for (FastaEntryInterface entry : list) {
			writer.write(entry);
		}
		writer.close();
	}
}
