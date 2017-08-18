package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.awt.Dimension;
import java.io.File;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;

public class MzmlStructureCharterTest {

	public static void main(String[] args) {
		File mzMLFile=new File("/Users/searleb/Documents/school/projects/pecandata/121115_BCS_HeLa_24mz_400_1000.mzML");
		mzMLFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/test/fast/121115_BCS_HeLa_24mz_400_1000.mzML");
		mzMLFile=new File("/Volumes/DataBackup/david_hawke/040117HeLa-dig-2ug-longDIA.mzML");
		// mzMLFile=new
		// File("/Volumes/DataBackup/david_hawke/040417HeLa-dig-2ug-longDIA.mzML");
		// mzMLFile=new
		// File("/Volumes/DataBackup/david_hawke/040117HeLa-dig-2ug-longerDIA.mzML");
		// mzMLFile=new
		// File("/Volumes/DataBackup/david_hawke/040417HeLa-dig-2ug-longDIA_2.mzML");
		// mzMLFile=new
		// File("/Volumes/DataBackup/david_hawke/040317HeLa-dig-2ug-longDIAthermo.mzML");
		mzMLFile=new File("/Users/searleb/Downloads/2017_January_23_envtstress_geoduck3.mzML");

		ChartPanel panel=MzmlStructureCharter.getStructureChart(mzMLFile);
		Charter.launchComponent(panel, mzMLFile.getName(), new Dimension(792, 612));
	}
}
