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
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class DilutionCurveFitterExample {

	
	public static void main(String[] args) throws Exception {
		final File outputDirectory=new File("/Users/searle.30/Downloads/ariana/Babies_first_BSA_curvefitting/");
		final File targetDirectory=new File(outputDirectory, "target");
		outputDirectory.mkdirs();
		targetDirectory.mkdirs();
		
		File dataFile=new File("/Users/searle.30/Downloads/ariana/Babies_first_BSA_quant.csv");
		File sampleOrganizationFile=new File("/Users/searle.30/Downloads/ariana/prm_sample_organization_BSA.csv");

		
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
}
