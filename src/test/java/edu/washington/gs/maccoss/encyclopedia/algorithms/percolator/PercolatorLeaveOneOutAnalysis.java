package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.io.LineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.LineParserMuscle;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class PercolatorLeaveOneOutAnalysis {
	private static final String DELIM = "\t";
	
	public static void main(String[] args) throws Exception {
		File testingPinFile=new File("/Users/searleb/Documents/encyclopedia/bugs/broken_v2_percolator/13Oct2020-Lumos_dlp_TPAD_CSF_DIA_chrlib_4mz_900to1000_ALLpoolLib_052.dia.features.txt");
		File fastaFile=new File("/Users/searleb/Documents/encyclopedia/bugs/broken_v2_percolator/uniprot_human_25apr2019.fasta");
		
		final String[] headers = new String[] { "xTandem", "xCorrLib", "xCorrModel", "dotProduct", "contrastAngle",
				"logit", "sumOfSquaredErrors", "numberOfMatchingPeaks", "numberOfMatchingPeaksAboveThreshold",
				"averageFragmentDeltaMasses", "isotopeDotProduct", "averageParentDeltaMass", "percentBlankOverMono",
				"numberPrecursorMatch", "sp", "maxLadderLength", "primary", "secondary", "evalue",
				"correlationToGaussian", "correlationToPrecursor", "isIntegratedSignal", "isIntegratedPrecursor",
				"numPeaksWithGoodCorrelation", "deltaRT", "ms1MassError", "ms2MassError", "numMissedCleavage",
				"precursorMz", "precursorMass", "RTinMin"};

		for (String leaveColumnOutName : headers) {
			BufferedReader in=new BufferedReader(new FileReader(testingPinFile));
			String header=in.readLine();
			String[] columnNames=header.split(DELIM, -1);

			int indexToLeaveOutTemp=-1;
			for (int i = 0; i < columnNames.length; i++) {
				if (leaveColumnOutName.equals(columnNames[i])) {
					indexToLeaveOutTemp=i;
					break;
				}
			}
			final int indexToLeaveOut=indexToLeaveOutTemp;

			System.out.println("============================ REMOVING "+leaveColumnOutName+" ===============================");
			in.close();
			
			//File tempPin=File.createTempFile("Percolator_", ".pin");
			File tempPin=new File(testingPinFile.getParentFile(), "remove_"+leaveColumnOutName+".pin");
			tempPin.deleteOnExit();

			FileWriter fileStream=new FileWriter(tempPin);
			PrintWriter out=new PrintWriter(fileStream);
			
			LineParserMuscle muscle=new LineParserMuscle() {
				
				@Override
				public void processRow(String row) {
					String[] data=row.split(DELIM, -1);
					
					for (int i = 0; i < data.length; i++) {
						if (i!=indexToLeaveOut) {
							if (i>0) out.print(DELIM);
							out.print(data[i]);
						}
					}
					out.println();
				}
				
				@Override
				public void cleanup() {
					
				}
			};
			LineParser.parseFile(testingPinFile, muscle);
			
			out.close();

			PercolatorExecutionData percolatorFiles=PercolatorExecutorTest.getPercolatorFiles(tempPin, fastaFile, SearchParameterParser.getDefaultParametersObject());

			final AminoAcidConstants aaConstants = new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());
			final float threshold = 0.01f;

			Pair<ArrayList<PercolatorPeptide>, Float> origpair=PercolatorExecutor.executePercolatorTSV(PercolatorVersion.v3p01, percolatorFiles, threshold, aaConstants, 2);
			
			System.out.println("----------------- FINAL "+leaveColumnOutName+" --> "+origpair.x.size()+" ("+origpair.y+") ------------------------------");
		}
	}

}
