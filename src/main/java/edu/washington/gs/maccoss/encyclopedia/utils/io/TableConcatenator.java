package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class TableConcatenator {
	private static final String DELIM = "\t";

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
	
	/**
	 * Parse TSV.
	 * Slower, but doesn't assume that tables all have the same headers (assumes headers from the first file). Also, only keeps best peptide based on a primary score
	 * Also assumes primary scores increase as they get better
	 * @param tables
	 */
	public static void concatenatePINTables(ArrayList<File> tables, File output, String primaryScore) throws IOException {
		String[] columnNames=null;
		HashMap<String, ScoredRow> dataset=new HashMap<>();

		for (File file : tables) {
			if (columnNames==null) {
				BufferedReader in=new BufferedReader(new FileReader(file));
				String header=in.readLine();
				columnNames=header.split(DELIM, -1);
				in.close();
			}
			
			TableParserMuscle muscle=new PINMuscle(primaryScore, dataset);
			TableParser.parseTSV(file, muscle);
		}
		
		ArrayList<ScoredRow> rows=new ArrayList<>(dataset.values());
		Collections.sort(rows);
		
		FileWriter fileStream=new FileWriter(output);
		PrintWriter out=new PrintWriter(fileStream);
		out.println(General.toString(columnNames, DELIM));
		for (ScoredRow row : rows) {
			Map<String, String> data=row.getRow();
			for (int i = 0; i < columnNames.length; i++) {
				if (i>0) out.print(DELIM);
				out.print(data.get(columnNames[i]));
			}
			out.println();
		}
		
		out.close();
	}
	
	public static class PINMuscle implements TableParserMuscle {
		RuntimeException error=null;
		private final String primaryScore;
		private final HashMap<String, ScoredRow> dataset; // only works because tableparser threading still runs just a single processing job
		
		
		public PINMuscle(String primaryScore, HashMap<String, ScoredRow> dataset) {
			super();
			this.primaryScore=primaryScore;
			this.dataset=dataset;
		}

		@Override
		public void processRow(Map<String, String> row) {
			String sequence=row.get("sequence");
			String primaryScoreValue=row.get(primaryScore);
			
			if (sequence==null) {
				error=new EncyclopediaException("Couldn't find sequence in PIN file!");
				throw error;
			} else if (primaryScoreValue==null) {
				error=new EncyclopediaException("Couldn't find sequence in primary score ("+primaryScoreValue+") file!");
				throw error; 
			}
			
			float primary=Float.parseFloat(primaryScoreValue);
			
			ScoredRow previous=dataset.get(sequence);
			if (previous==null) {
				dataset.put(sequence, new ScoredRow(sequence, primary, row));
			} else {
				if (previous.getScore()<primary) {
					dataset.put(sequence, new ScoredRow(sequence, primary, row));
				}
			}
		}
		
		@Override
		public void cleanup() {
		}
		
		public Exception getError() {
			return error;
		}
	}
	
	public static class ScoredRow implements Comparable<ScoredRow> {
		private final String peptide;
		private final float score;
		private final Map<String, String> row;
		
		public ScoredRow(String peptide, float score, Map<String, String> row) {
			this.peptide=peptide;
			this.score=score;
			this.row=row;
		}

		@Override
		public int compareTo(ScoredRow o) {
			if (o==null) return 1;
			int c=peptide.compareTo(o.peptide);
			if (c!=0) return c;
			
			c=Float.compare(score, o.score);
			return c;
		}
		
		@Override
		public int hashCode() {
			return Float.hashCode(score);
		}
		
		public float getScore() {
			return score;
		}
		
		public Map<String, String> getRow() {
			return row;
		}
	}
}
