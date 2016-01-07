package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class TableConcatenator {
	/**
	 * Assumes that all tables have the same number and location of columns.
	 * Literally just does a direct copy of everything after the header.
	 * 
	 * @param tables
	 */
	public static void concatenateTables(ArrayList<File> tables, File output) throws IOException {
		FileWriter fileStream=new FileWriter(output);
		BufferedWriter out=new BufferedWriter(fileStream);

		boolean isFirst=true;
		for (File file : tables) {
			BufferedReader in=new BufferedReader(new FileReader(file));
			if (!isFirst) {
				// drop first line
				in.readLine();
			}

			char[] buff=new char[1024];

			int n;
			while (-1!=(n=in.read(buff))) {
				out.write(buff, 0, n);
			}
			in.close();
			isFirst=false;
		}
		
		out.close();
	}
}
