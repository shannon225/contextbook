package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.StringUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.io.LineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.LineParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.math.FloatPair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearDiscriminantAnalysis;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearInterpolatedFunction;
import edu.washington.gs.maccoss.encyclopedia.utils.math.PivotTableGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RunningMedianWarper;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.math.randomforest.RocPlot;
import gnu.trove.list.array.TFloatArrayList;

public class MProphetFeatureReader {
	private static final String DELIM = "\t";

	public static MProphetDataset parseFeatureFile(File file, MProphetExecutionData settings) throws EncyclopediaException {
		String[] columnNames=null;
		int idIndex=0;
		int labelIndex=0;
		int scanNrIndex=0;
		int sequenceIndex=0;
		int proteinIndex=0;
		
		boolean[] isFeature=null;
		ArrayList<String> featureNames=new ArrayList<>();

		Logger.logLine("Parsing header for input file "+file.getName());
		try {
			BufferedReader in=new BufferedReader(new FileReader(file));
			String header=in.readLine();
			columnNames=header.split(DELIM, -1);
			isFeature=new boolean[columnNames.length];
			for (int i = 0; i < columnNames.length; i++) {
				if ("id".equals(columnNames[i])) {
					idIndex=i;
				} else if ("Label".equals(columnNames[i])) {
					labelIndex=i; // either TD or Label
				} else if ("TD".equals(columnNames[i])) {
					labelIndex=i; // either TD or Label
				} else if ("ScanNr".equals(columnNames[i])) {
					scanNrIndex=i;
				} else if ("sequence".equals(columnNames[i])) {
					sequenceIndex=i;
				} else if ("Proteins".equals(columnNames[i])) {
					proteinIndex=i;
				} else if ("protein".equals(columnNames[i])) {
					proteinIndex=i;
				} else if ("pepLength".equals(columnNames[i])) {
					// skip
				} else if ("charge1".equals(columnNames[i])) {
					// skip
				} else if ("charge2".equals(columnNames[i])) {
					// skip
				} else if ("charge3".equals(columnNames[i])) {
					// skip
				} else if ("charge4".equals(columnNames[i])) {
					// skip
				} else if ("precursorMass".equals(columnNames[i])) {
					// skip
				} else if ("RTinMin".equals(columnNames[i])) {
					// skip
				} else if ("midTime".equals(columnNames[i])) {
					// skip
				} else if ("numberOfMatchingPeaksAboveThreshold".equals(columnNames[i])) {
					// skip
				} else if ("primary".equals(columnNames[i])) {
					// Primary score can get overloaded, so drop!
				} else {
					featureNames.add(columnNames[i]);
					isFeature[i]=true;
				}
			}
			
			in.close();
			
			String sortingScoreString=null;
			
			// prefer last in list (xCorrModel)
			if (StringUtils.contains(columnNames, "peakZScore")) sortingScoreString="peakZScore";  
			if (StringUtils.contains(columnNames, "peakBGScore")) sortingScoreString="peakBGScore"; 
			if (StringUtils.contains(columnNames, "peakBGSScore")) sortingScoreString="peakBGSScore";  
			if (StringUtils.contains(columnNames, "xTandem")) sortingScoreString="xTandem";  
			if (StringUtils.contains(columnNames, "scribe")) sortingScoreString="scribe";  
			if (StringUtils.contains(columnNames, "xCorrLib")) sortingScoreString="xCorrLib";  
			if (StringUtils.contains(columnNames, "xCorrModel")) sortingScoreString="xCorrModel";  
			//if (StringUtils.contains(columnNames, "primary")) sortingScoreString="primary";  

			if (sortingScoreString==null) {
				throw new EncyclopediaException("Can't parse sorting score from header from ["+header+"]");
			} else {
				Logger.logLine("Using ["+sortingScoreString+"] as the initial score.");
			}
			int sortingScoreIndex=0;
			for (int i = 0; i < columnNames.length; i++) {
				if (columnNames[i].equals(sortingScoreString)) {
					sortingScoreIndex=i;
					break;
				}
			}
			Logger.logLine("Found indicies for "+featureNames.size()+" features: ["+General.toString(featureNames)+"]");

			final int idIndexFinal=idIndex;
			final int labelIndexFinal=labelIndex;
			final int sequenceIndexFinal=sequenceIndex;
			final int proteinIndexFinal=proteinIndex;
			final boolean[] isFeatureFinal=isFeature;
			final String[] columnNamesFinal=columnNames;
			
			ArrayList<MProphetData> peptideData=new ArrayList<>();
			LineParserMuscle muscle = new LineParserMuscle() {
				boolean isFirst=true;
				
				@Override
				public void processRow(String row) {
					if (isFirst) {
						// skip header
						isFirst=false;
						return;
					}
					String[] values=row.split(DELIM, -1);
					boolean isDecoy;
					
					try {
						isDecoy=Integer.parseInt(values[labelIndexFinal])<0;
					} catch (Exception e) {
						Logger.errorLine("Error parsing ["+values[labelIndexFinal]+"] as "+columnNamesFinal[labelIndexFinal]+" (index "+labelIndexFinal+")!");
						Logger.errorException(e);
						throw e;
					}
					
					TFloatArrayList features=new TFloatArrayList();
					for (int j = 0; j < isFeatureFinal.length; j++) {
						if (isFeatureFinal[j]) {
							try {
								features.add(Float.parseFloat(values[j]));
							} catch (Exception e) {
								Logger.errorLine("Error parsing ["+values[j]+"] as "+columnNamesFinal[j]+" (index "+j+")!");
								Logger.errorException(e);
								throw e;
							}
						}
					}
					String peptideModSeq=PercolatorPeptide.getPeptideSequence(values[idIndexFinal]);
					peptideData.add(new MProphetData(values[idIndexFinal], peptideModSeq, values[proteinIndexFinal], features.toArray(), isDecoy));
					
				}
				
				@Override
				public void cleanup() {
				}
			};
			LineParser.parseFile(settings.getInputTSV(), muscle);

			return new MProphetDataset(featureNames, sortingScoreIndex, peptideData);
			
		} catch (Throwable t) {
			Logger.errorLine("Error performing mProphet!");
			Logger.errorException(t);

			throw new EncyclopediaException(t);
		}
	}
}
