package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.awt.Dimension;
import java.awt.GridLayout;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.JPanel;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;

public class ChromatogramExtractorTest {
	public static void main(String args[]) throws Exception {
		File diaFile=new File("/Users/searleb/Documents/phospho_localization/data/hela/110515_bcs_hela_phospho_starved_20mz_500_900.dia");

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		HashMap<String, FragmentationModel> availableModels=new HashMap<String, FragmentationModel>();
		byte charge=(byte)2;
		
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
		String peptideModSeq1="KT[+79.966331]APTLSPEHWK";
		FragmentationModel model1=PeptideUtils.getPeptideModel(peptideModSeq1, parameters.getAAConstants());
		FragmentIon[] ionTypes1=model1.getPrimaryIonObjects(FragmentationType.CID, charge, false);
		availableModels.put(peptideModSeq1, model1);

		String peptideModSeq2="KTAPTLS[+79.966331]PEHWK";
		FragmentationModel model2=PeptideUtils.getPeptideModel(peptideModSeq2, parameters.getAAConstants());
		FragmentIon[] ionTypes2=model2.getPrimaryIonObjects(FragmentationType.CID, charge, false);
		availableModels.put(peptideModSeq2, model2);

		FragmentIon[] unique1=PhosphoLocalizer.getUniqueFragmentIons(peptideModSeq1, charge, availableModels, parameters);
		FragmentIon[] unique2=PhosphoLocalizer.getUniqueFragmentIons(peptideModSeq2, charge, availableModels, parameters);
		
		ArrayList<FragmentScan> stripes=stripefile.getStripes(737.858763, 62f*60, 67f*60, false);

		ArrayList<Spectrum> spectra=new ArrayList<Spectrum>();
		for (FragmentScan stripe : stripes) {
			spectra.add(stripe);
		}

		JPanel panel=new JPanel(new GridLayout(0, 2));

		ArrayList<ChartPanel> charts=new ArrayList<ChartPanel>();
		double globalMaxY=0.0;
		double[] ppms=new double[] {5, 10, 16.67, 20, 50};
		for (int i=0; i<ppms.length; i++) {
			boolean displayLegend=i==ppms.length-1;
			MassTolerance tol=new MassTolerance(ppms[i]);
			
			HashMap<FragmentIon, XYTrace> map1=ChromatogramExtractor.extractFragmentChromatograms(tol, ionTypes1, spectra, null, GraphType.dashedline);
			HashMap<FragmentIon, XYTrace> uniqueMap1=ChromatogramExtractor.extractFragmentChromatograms(tol, unique1, spectra, null, GraphType.boldline);
			map1.putAll(uniqueMap1);
			ChartPanel chart1=Charter.getChart(peptideModSeq1+" Retention Time ("+tol.toString()+")", "Intensity", displayLegend, map1.values().toArray(new XYTrace[map1.size()]));
			panel.add(chart1);
			charts.add(chart1);
			
			HashMap<FragmentIon, XYTrace> map2=ChromatogramExtractor.extractFragmentChromatograms(tol, ionTypes2, spectra, null, GraphType.dashedline);
			HashMap<FragmentIon, XYTrace> uniqueMap2=ChromatogramExtractor.extractFragmentChromatograms(tol, unique2, spectra, null, GraphType.boldline);
			map2.putAll(uniqueMap2);
			ChartPanel chart2=Charter.getChart(peptideModSeq2+" Retention Time ("+tol.toString()+")", "Intensity", displayLegend, map2.values().toArray(new XYTrace[map2.size()]));
			panel.add(chart2);
			charts.add(chart2);

			for (XYTrace trace : uniqueMap1.values()) {
				globalMaxY=Math.max(globalMaxY, trace.getMaxY());
			}
			for (XYTrace trace : uniqueMap2.values()) {
				globalMaxY=Math.max(globalMaxY, trace.getMaxY());
			}
		}
		
		if (globalMaxY>0.0) {
			globalMaxY=globalMaxY*1.05;
			for (ChartPanel chartPanel : charts) {
				chartPanel.getChart().getXYPlot().getRangeAxis().setUpperBound(globalMaxY);
			}
		}

		Charter.launchComponent(panel, "KTAPTLSPEHWK [+80.0]", new Dimension(1900, 1030));
		
	}
}
