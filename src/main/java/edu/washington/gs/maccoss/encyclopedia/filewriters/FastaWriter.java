package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;

import org.apache.commons.lang3.text.WordUtils;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class FastaWriter {
	public static void writeFasta(File f, ArrayList<FastaEntryInterface> list) {
		PrintWriter writer=null;
		try {
			writer=new PrintWriter(f, "UTF-8");
			for (FastaEntryInterface entry : list) {
				writer.print('>');
				writer.println(entry.getAccession());
				writer.println(WordUtils.wrap(entry.getSequence(), 80));
			}
		} catch (FileNotFoundException e) {
			Logger.logException(e);
		} catch (UnsupportedEncodingException e) {
			Logger.logException(e);
		} finally {
			if (writer!=null) {
				writer.close();
			}
		}
	}
}
