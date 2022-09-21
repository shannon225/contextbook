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
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
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
	 * Only keeps best peptide based on a primary score. Assumes primary scores increase as they get better!
	 * @param tables
	 */
	public static void concatenatePINTables(ArrayList<File> tables, File output, String primaryScore) throws IOException {
		String[] columnNames=null;
		int sequenceIndex=0;
		int primaryScoreIndex=0;
		
		HashMap<String, ScoredRow> dataset=new HashMap<>();

		for (File file : tables) {
			Logger.logLine("Parsing Percolator input file "+file.getName());
			
			if (columnNames==null) {
				BufferedReader in=new BufferedReader(new FileReader(file));
				String header=in.readLine();
				columnNames=header.split(DELIM, -1);
				for (int i = 0; i < columnNames.length; i++) {
					if ("sequence".equals(columnNames[i])) {
						sequenceIndex=i;
					} else if (primaryScore.equals(columnNames[i])) {
						primaryScoreIndex=i;
					}
				}
				in.close();
				Logger.logLine("Found indicies for sequence: "+sequenceIndex+" and "+primaryScore+": "+primaryScoreIndex);
			}
			
			LineParserMuscle muscle=new PINMuscle(dataset, sequenceIndex, primaryScoreIndex);
			LineParser.parseFile(file, muscle);
		}
		
		ArrayList<ScoredRow> rows=new ArrayList<>(dataset.values());
		Collections.sort(rows);
		Logger.logLine("Found "+rows.size()+" total peptides, writing to new Percolator input file...");
		
		FileWriter fileStream=new FileWriter(output);
		PrintWriter out=new PrintWriter(fileStream);
		out.println(General.toString(columnNames, DELIM));
		for (ScoredRow row : rows) {
			out.println(row.getRow());
		}
		
		out.close();
	}
	
	public static class PINMuscle implements LineParserMuscle {
		RuntimeException error=null;
		private final HashMap<String, ScoredRow> dataset; // only works because tableparser threading still runs just a single processing job
		private final int sequenceIndex;
		private final int primaryScoreIndex;
		private boolean isFirst=true;
		
		public PINMuscle(HashMap<String, ScoredRow> dataset, int sequenceIndex, int primaryScoreIndex) {
			super();
			this.dataset=dataset;
			this.sequenceIndex=sequenceIndex;
			this.primaryScoreIndex=primaryScoreIndex;
		}

		@Override
		public void processRow(String row) {
			if (isFirst) {
				// skip header
				isFirst=false;
				return;
			}
			String[] data=row.split(DELIM, -1);
			
			String sequence=data[sequenceIndex];
			String primaryScoreValue=data[primaryScoreIndex];
			
			if (sequence==null) {
				error=new EncyclopediaException("Couldn't find sequence in PIN file!");
				throw error;
			} else if (primaryScoreValue==null) {
				error=new EncyclopediaException("Couldn't find primary score ("+primaryScoreValue+")! Index: "+primaryScoreIndex+", Row: "+row);
				throw error; 
			}
			
			float primary;
			try {
				primary=Float.parseFloat(primaryScoreValue);
			} catch (NumberFormatException nfe) {
				error=new EncyclopediaException("Couldn't parse primary score ("+primaryScoreValue+")! Index: "+primaryScoreIndex+", Row: "+row);
				throw error; 
			}
			
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
		private final String row;
		
		public ScoredRow(String peptide, float score, String row) {
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
		
		public String getRow() {
			return row;
		}
	}
}
