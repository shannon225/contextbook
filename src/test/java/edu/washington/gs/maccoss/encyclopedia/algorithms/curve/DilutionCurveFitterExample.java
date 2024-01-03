package edu.washington.gs.maccoss.encyclopedia.algorithms.curve;

import java.awt.Dimension;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.algorithms.curve.DilutionCurveFitter.FitPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class DilutionCurveFitterExample {

	
	public static void main2(String[] args) throws Exception {
		final File outputDirectory=new File("/Users/searleb/Downloads/LOD_LOQ_yi_sampleprep_titration_updated/curve_fitting");
		final File targetDirectory=new File(outputDirectory, "target");
		outputDirectory.mkdirs();
		targetDirectory.mkdirs();
		
		File dataFile=new File("/Users/searleb/Downloads/LOD_LOQ_yi_sampleprep_titration_updated/2022_09_26_cell_titration_exp_quant.elib.peptides.txt");
		File sampleOrganizationFile=new File("/Users/searleb/Downloads/LOD_LOQ_yi_sampleprep_titration_updated/prm_sample_organization_tcell_sample_prep.csv");

		
		Pair<ArrayList<ScoredObject<String>>, Map<String, TObjectFloatHashMap<String>>> concentrationPair=DilutionCurveFitter.getExpectedConcentrationsFromCSV(sampleOrganizationFile);
		final ArrayList<ScoredObject<String>> expectedConcentrations=concentrationPair.x;
		final Map<String, TObjectFloatHashMap<String>> unknowns=concentrationPair.y;
		
		for (ScoredObject<String> scoredObject : expectedConcentrations) {
			System.out.println("Expected: "+scoredObject);
		}
		for (Entry<String, TObjectFloatHashMap<String>> entry : unknowns.entrySet()) {
			System.out.println("Unknown: "+entry.getKey()+"\t"+General.toString(entry.getValue().keys()));
		}
		
		final float[] expected = DilutionCurveFitter.adjustForZeroConcentrations(expectedConcentrations);

		final ArrayList<FitPeptide> fitPeptides=DilutionCurveFitter.fitCurves(outputDirectory, dataFile, expectedConcentrations, expected, null, true);
		final HashMap<String, Map<String, TObjectFloatHashMap<String>>> unknownData=DilutionCurveFitter.extractUnknowns(dataFile, unknowns, null);
		
		for (FitPeptide fit : fitPeptides) {
			System.out.println("Peptide: "+fit.getPeptideModSeq());
			
			Map<String, TObjectFloatHashMap<String>> data=unknownData.get(fit.getPeptideModSeq());
			
			ChartPanel panel=DilutionCurveFitter.graph(fit.getPeptideModSeq(), fit.getExpectedRelativeIntensities(), fit.getActualRelativeIntensities(), fit.getBestFit(), Optional.ofNullable(data));
			Charter.writeAsPDF(panel.getChart(), new File(targetDirectory, fit.getPeptideModSeq()+".pdf"), new Dimension(400, 300));
		}
	}
	
	public static void main(String[] args) throws Exception {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		File inputDirectory=new File("/Users/searleb/Documents/students/ariana/");
		File outputDirectory=new File("/Users/searleb/Documents/students/ariana/curvefitting_combined_30k_testing_10p/");
		File dataFile=new File(inputDirectory, "dimethyl_combined_library_titration_curve.elib.peptides.txt");
		File sampleOrganizationFile=new File(inputDirectory, "exploris_21pt_curve_samples_list.csv");
		File libraryFile=new File(inputDirectory, "v3_library_exploris.elib");
		File rtAlignFile=new File(inputDirectory, "2024_01_02_mouse_WT_CD8Tcell_250ng_16mzst_DIA_auorora_01.mzML.elib");

		for (boolean isDeepAssay : new boolean [] {true}) {
			AbstractDilutionCurveFittingParameters fittingParams;
			if (isDeepAssay) {
				fittingParams=new DilutionCurveParameters(true, true);
			} else {
				fittingParams=new DilutionCurveParameters(true, false);
			}
			
			DilutionCurveFitter.generateAssayFromCurves(params, outputDirectory, dataFile, sampleOrganizationFile, libraryFile, rtAlignFile, fittingParams);
		}
	}
	
	public static String[] targets = new String[] { "A2APV2", "B2RQC6", "E9PVX6", "E9Q394", "O08538", "O08573",
			"O08738", "O35188", "O35205", "O35228", "O35405", "O35608", "O35714", "O35892", "O54707", "O54824",
			"O54839", "O54907", "O55038", "O55237", "O70460", "O88324", "O88713", "O89093", "O89110", "P01132",
			"P01575", "P01580", "P01582", "P01731", "P01831", "P02802", "P03958", "P04187", "P04202", "P04351",
			"P05064", "P06332", "P06804", "P07141", "P07750", "P08113", "P08505", "P08882", "P08884", "P09793",
			"P10148", "P10417", "P10518", "P10820", "P10852", "P10855", "P11032", "P11440", "P12850", "P13379",
			"P13864", "P14097", "P14901", "P15247", "P15379", "P15655", "P15702", "P16045", "P17515", "P18337",
			"P18340", "P18572", "P20029", "P20334", "P26262", "P27512", "P27548", "P27814", "P28352", "P28654",
			"P28667", "P29391", "P31041", "P31240", "P34960", "P35550", "P35918", "P37172", "P37217", "P40224",
			"P41047", "P41272", "P42227", "P43432", "P46978", "P47741", "P48346", "P49764", "P49945", "P50228",
			"P50592", "P52293", "P54227", "P55194", "P55821", "P56203", "P56528", "P60324", "P63089", "P70313",
			"P70677", "P83555", "P86176", "Q00417", "Q00731", "Q01320", "Q02242", "Q02858", "Q03366", "Q07763",
			"Q08048", "Q10738", "Q149L7", "Q3TDQ1", "Q3U7R1", "Q60715", "Q60837", "Q60865", "Q60875", "Q61176",
			"Q61790", "Q61823", "Q62293", "Q62351", "Q64261", "Q66JW3", "Q6A0A9", "Q6NZJ6", "Q7TMM9", "Q7TNG5",
			"Q7TST5", "Q80WS3", "Q80XI3", "Q8BRF7", "Q8BVZ5", "Q8C147", "Q8C156", "Q8C196", "Q8C522", "Q8C567",
			"Q8JZK9", "Q8K209", "Q8VCW8", "Q8VHB5", "Q8VIM0", "Q91V12", "Q922F4", "Q99J36", "Q99JB6", "Q99JR5",
			"Q99KK2", "Q99PA5", "Q9CQE5", "Q9CQU0", "Q9CR75", "Q9D659", "Q9D6F9", "Q9DC29", "Q9EP73", "Q9EPE9",
			"Q9EPU5", "Q9ERD8", "Q9ES17", "Q9ET39", "Q9JHH5", "Q9JHJ8", "Q9JL26", "Q9JLV6", "Q9QYH9", "Q9R1Q7",
			"Q9WTK5", "Q9WUL5", "Q9WVS0", "Q9Z0D9", "Q9Z121", };

	public static class DilutionCurveParameters implements AbstractDilutionCurveFittingParameters {
		final int numberOfRTAnchors=0;
		final int maxNumberPeptidesPerProtein;
		final int targetTotalNumberOfPeptides=300; // remember to subtract off anchors (total is 160 peptides)
		final float windowInMin=5f; // in minutes!
		final float minCVForAnchors=0.05f;
		final float minCVForBadAnchors=0.75f;
		final int assayMaxDensity;
		final boolean requireAlignmentRT=false; // turn off for fitting against PRM
		final boolean useLineNoise=true; // newer versions should set this to "true"
		final boolean untargeted;
		
		public DilutionCurveParameters(boolean untargeted, boolean isdeep) {
			super();
			this.untargeted = untargeted;
			if (isdeep) {
				maxNumberPeptidesPerProtein=3;
				assayMaxDensity=10;
			} else {
				maxNumberPeptidesPerProtein=5;
				assayMaxDensity=20;
			}
		}
		
		public float getWindowInMin() {
			return windowInMin;
		}
		public float getWindowInMin(float rtInSec) {
			return windowInMin;
		}

		public int getNumberOfRTAnchors() {
			return numberOfRTAnchors;
		}

		public int getMaxNumberPeptidesPerProtein() {
			return maxNumberPeptidesPerProtein;
		}

		public int getTargetTotalNumberOfPeptides() {
			return targetTotalNumberOfPeptides;
		}

		public float getMinCVForAnchors() {
			return minCVForAnchors;
		}

		public float getMinCVForBadAnchors() {
			return minCVForBadAnchors;
		}

		public int getAssayMaxDensity() {
			return assayMaxDensity;
		}

		@Override
		public boolean isTargetedProtein(String accession) {
			if (untargeted) return true;
			
			for (int i = 0; i < targets.length; i++) {
				if (accession.indexOf(targets[i])>=0) return true;
			}
			return false;
		}
		
		@Override
		public boolean isEliminatedPeptide(String peptideModSeq) {
			return false;
		}

		public boolean isRequireAlignmentRT() {
			return requireAlignmentRT;
		}

		public boolean isUseLineNoise() {
			return useLineNoise;
		}
	}
}
