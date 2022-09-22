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
import java.util.function.BiConsumer;

import edu.washington.gs.maccoss.encyclopedia.filewriters.ScoringResultsToTSVConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedianDouble;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.KDE;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class TableConcatenator {
	private static final float TARGET_CONSISTENCY_PERCENTAGE = 0.5f;
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
		if (tables.size()==0) return;
		
		String[] columnNames=null;
		int sequenceIndex=0;
		int primaryScoreIndex=0;
		int deltaRTIndex=-1;
		int midTimeIndex=-1;
		
		// get column names from first file
		BufferedReader in=new BufferedReader(new FileReader(tables.get(0)));
		String header=in.readLine();
		columnNames=header.split(DELIM, -1);
		String deltaRTName = ScoringResultsToTSVConsumer.deltaRTName;
		for (int i = 0; i < columnNames.length; i++) {
			if ("sequence".equals(columnNames[i])) {
				sequenceIndex=i;
			} else if (deltaRTName.equals(columnNames[i])) {
				deltaRTIndex=i;
			} else if ("midTime".equals(columnNames[i])) {
				midTimeIndex=i; // only Pecan uses this, and it doesn't have deltaRT
			} else if (primaryScore.equals(columnNames[i])) {
				primaryScoreIndex=i;
			}
		}
		in.close();
		if (deltaRTIndex==-1&&midTimeIndex>=0) {
			deltaRTIndex=midTimeIndex;
			deltaRTName="midTime";
		}
		Logger.logLine("Found indicies for sequence: "+sequenceIndex+", "+primaryScore+": "+primaryScoreIndex+", and "+deltaRTName+": "+deltaRTIndex);

		// calculate deltaRT distributions
		HashMap<String, TFloatArrayList> scoresBySequence=new HashMap<>();
		HashMap<String, TFloatArrayList> deltaRTsBySequence=new HashMap<>();
		final HashMap<String, DeltaRTData> deltaRTDataBySequence=new HashMap<>();
		boolean scoreRT = deltaRTIndex>=0;
		if (scoreRT) {
			for (File file : tables) {
				Logger.logLine("Getting delta retention times for input file "+file.getName());
				DeltaRTMuscle muscle=new DeltaRTMuscle(deltaRTsBySequence, scoresBySequence, sequenceIndex, primaryScoreIndex, deltaRTIndex);
				LineParser.parseFile(file, muscle);
			}
			deltaRTsBySequence.forEach(new BiConsumer<String, TFloatArrayList>() {
				@Override
				public void accept(String t, TFloatArrayList u) {
					deltaRTDataBySequence.put(t, new DeltaRTData(u, scoresBySequence.get(t)));
				}
			});
		}
		
		// parse files with muscle to select best row
		HashMap<String, ScoredRow> dataset=new HashMap<>();
		for (File file : tables) {
			Logger.logLine("Parsing Percolator input file "+file.getName());
			LineParserMuscle muscle=new PINMuscle(dataset, sequenceIndex, primaryScoreIndex, deltaRTIndex, deltaRTDataBySequence);
			LineParser.parseFile(file, muscle);
		}
		
		// collapse final list
		ArrayList<ScoredRow> rows=new ArrayList<>(dataset.values());
		Collections.sort(rows);
		Logger.logLine("Found "+rows.size()+" total peptides, writing to new Percolator input file...");
		
		FileWriter fileStream=new FileWriter(output);
		PrintWriter out=new PrintWriter(fileStream);
		if (scoreRT) {
			columnNames=General.insert(columnNames, 3, "consistencyScore", "absDeltaDeltaRT", "irqDeltaRT");
		}
		out.println(General.toString(columnNames, DELIM));
		for (ScoredRow row : rows) {
			out.println(row.getRow(scoreRT));
		}
		
		out.close();
	}
	
	protected static class DeltaRTData {
		private final float modeDeltaRT;
		private final float irqDeltaRT;
		private final float top20PercentPrimaryScore;
		public DeltaRTData(TFloatArrayList deltaRTs, TFloatArrayList primaryScores) {
			KDE distribution=new KDE(General.toDoubleArray(deltaRTs.toArray()), 1.0);
			this.modeDeltaRT = (float)distribution.getMode();
			this.irqDeltaRT = QuickMedian.iqr(deltaRTs.toArray());
			
			// need at least 3 samples to be "consistent"
			float consistencyPercentage=Math.max(TARGET_CONSISTENCY_PERCENTAGE, 3.0f/primaryScores.size());
			this.top20PercentPrimaryScore=QuickMedian.select(primaryScores.toArray(), consistencyPercentage);
		}
	}
	
	protected static class DeltaRTMuscle implements LineParserMuscle {
		RuntimeException error=null;
		private final HashMap<String, TFloatArrayList> deltaRTList; // only works because tableparser threading still runs just a single processing job
		private final HashMap<String, TFloatArrayList> scoreList; // only works because tableparser threading still runs just a single processing job
		private final int sequenceIndex;
		private final int primaryScoreIndex;
		private final int deltaRTIndex;
		private boolean isFirst=true;
		public DeltaRTMuscle(HashMap<String, TFloatArrayList> deltaRTList, HashMap<String, TFloatArrayList> scoreList, int sequenceIndex, int primaryScoreIndex, int deltaRTIndex) {
			super();
			this.deltaRTList=deltaRTList;
			this.scoreList=scoreList;
			this.sequenceIndex=sequenceIndex;
			this.primaryScoreIndex=primaryScoreIndex;
			this.deltaRTIndex=deltaRTIndex;
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
			String deltaRTString=data[deltaRTIndex];
			
			if (sequence==null) {
				error=new EncyclopediaException("Couldn't find sequence in PIN file!");
				throw error;
			} else if (deltaRTString==null) {
				error=new EncyclopediaException("Couldn't find deltaRT ("+deltaRTString+")! Index: "+deltaRTIndex+", Row: "+row);
				throw error; 
			}
			
			float primary;
			try {
				primary=Float.parseFloat(primaryScoreValue);
			} catch (NumberFormatException nfe) {
				error=new EncyclopediaException("Couldn't parse primary score ("+primaryScoreValue+")! Index: "+primaryScoreIndex+", Row: "+row);
				throw error; 
			}
			
			float deltaRT;
			try {
				deltaRT=Float.parseFloat(deltaRTString);
			} catch (NumberFormatException nfe) {
				error=new EncyclopediaException("Couldn't parse deltaRT ("+deltaRTString+")! Index: "+deltaRTIndex+", Row: "+row);
				throw error; 
			}
			
			TFloatArrayList list=scoreList.get(sequence);
			if (list==null) {
				list=new TFloatArrayList();
				scoreList.put(sequence, list);
			}
			list.add(primary);
			
			list=deltaRTList.get(sequence);
			if (list==null) {
				list=new TFloatArrayList();
				deltaRTList.put(sequence, list);
			}
			list.add(deltaRT);
		}
		
		@Override
		public void cleanup() {
		}
		
		public Exception getError() {
			return error;
		}
	}
	
	protected static class PINMuscle implements LineParserMuscle {
		RuntimeException error=null;
		private final HashMap<String, ScoredRow> dataset; // only works because tableparser threading still runs just a single processing job
		private final int sequenceIndex;
		private final int primaryScoreIndex;
		private final int deltaRTIndex;
		private final HashMap<String, DeltaRTData> deltaRTDataBySequence;
		private boolean isFirst=true;
		
		public PINMuscle(HashMap<String, ScoredRow> dataset, int sequenceIndex, int primaryScoreIndex, int deltaRTIndex, HashMap<String, DeltaRTData> deltaRTDataBySequence) {
			super();
			this.dataset=dataset;
			this.sequenceIndex=sequenceIndex;
			this.primaryScoreIndex=primaryScoreIndex;
			this.deltaRTIndex=deltaRTIndex;
			this.deltaRTDataBySequence=deltaRTDataBySequence;
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
			
			DeltaRTData deltaData=deltaRTDataBySequence.get(sequence);
			float deltaDeltaRT, irqDeltaRT, top20PercentPrimaryScore;
			if (deltaData!=null) {
				String deltaRTString=data[deltaRTIndex];
				
				float deltaRT;
				try {
					deltaRT=Float.parseFloat(deltaRTString);
				} catch (NumberFormatException nfe) {
					error=new EncyclopediaException("Couldn't parse deltaRT ("+deltaRTString+")! Index: "+deltaRTIndex+", Row: "+row);
					throw error; 
				}
				
				deltaDeltaRT=Math.abs(deltaData.modeDeltaRT-deltaRT);
				irqDeltaRT=deltaData.irqDeltaRT;
				top20PercentPrimaryScore=deltaData.top20PercentPrimaryScore;
			} else {
				deltaDeltaRT=0.0f;
				irqDeltaRT=0.0f;
				top20PercentPrimaryScore=0.0f;
			}
			
			ScoredRow previous=dataset.get(sequence);
			if (previous==null) {
				dataset.put(sequence, new ScoredRow(sequence, primary, primary, deltaDeltaRT, irqDeltaRT, top20PercentPrimaryScore, row));
			} else {
				// if the new score is either 10% better than the previous max score OR the delta is closer (and it's within 10% of the previous max score)
				if (previous.getTruemaxScore()<primary*0.9f||(previous.getTruemaxScore()*0.9f<primary&&previous.getDeltaDeltaRT()>deltaDeltaRT)) {
					dataset.put(sequence, new ScoredRow(sequence, primary, Math.max(primary, previous.getTruemaxScore()), deltaDeltaRT, irqDeltaRT, top20PercentPrimaryScore, row));
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
	
	protected static class ScoredRow implements Comparable<ScoredRow> {
		private final String peptide;
		private final float truemaxScore;
		private final float score;
		private final float deltaDeltaRT;
		private final float irqDeltaRT;
		private final float top20PercentPrimaryScore;
		private final String row;
		
		public ScoredRow(String peptide, float score, float truemaxScore, float deltaDeltaRT, float irqDeltaRT, float top20PercentPrimaryScore, String row) {
			this.peptide=peptide;
			this.score=score;
			this.truemaxScore=truemaxScore;
			this.deltaDeltaRT=deltaDeltaRT;
			this.irqDeltaRT=irqDeltaRT;
			this.top20PercentPrimaryScore=top20PercentPrimaryScore;
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
		
		public float getTruemaxScore() {
			return truemaxScore;
		}
		
		public float getDeltaDeltaRT() {
			return deltaDeltaRT;
		}

		public String getRow(boolean addAuxScores) {
			// row + top20PercentPrimaryScore, absDeltaDeltaRT, irqDeltaRT
			if (addAuxScores) {
				String[] values=row.split(DELIM, -1);
				String[] insert=General.insert(values, 3, Float.toString(top20PercentPrimaryScore), Float.toString(deltaDeltaRT), Float.toString(irqDeltaRT));
				return General.toString(insert, DELIM);
			} else {
				return row;
			}
		}
	}
}
