package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

import javax.swing.JFrame;
import javax.swing.JLabel;

import org.jzy3d.chart.ChartLauncher;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;

public class Spectra3dPanelTest {
	private static final int RT_MARGIN=0;

	public static void main(String[] args) throws Exception {
		SearchParameters params=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), DataAcquisitionType.OVERLAPPING_DIA, false, true);

		File file=new File("/Users/searleb/Documents/school/perspective/rawfiles/2017dec27_variable_dia_6b_rep1.mzML");
		StripeFileInterface raw=StripeFileGenerator.getFile(file, params, true);
		
		String peptideModSeq="VFSGLVSTGLK";
		byte charge=2;
		float targetRT=58.3f;
		
		peptideModSeq="APILIATDVASR";
		targetRT=58.4f;
		
		peptideModSeq="IMNVIGEPIDER";
		targetRT=58.55f;
		
		peptideModSeq="MNVLADALK";
		targetRT=58.95f;
		
		peptideModSeq="LVLVGDGGTGK";
		targetRT=44.15f;
		
		peptideModSeq="FADLSEAANR"; 
		targetRT=40.8f;
		
		peptideModSeq="LSGGLGAGSC[+57.0214635]R";
		targetRT=28.8f;
		
		
		float rtInSecStart=targetRT*60f-30;
		float rtInSecStop=targetRT*60f+15;
		
		FragmentationModel model=PeptideUtils.getPeptideModel(peptideModSeq, params.getAAConstants());
		FragmentIon[] ions=model.getPrimaryIonObjects(params.getFragType(), charge, true);
		double targetMz=model.getChargedMass(charge);
		
		
		ArrayList<Stripe> stripes=raw.getStripes(targetMz, rtInSecStart-RT_MARGIN, rtInSecStop+RT_MARGIN, false);
		Collections.sort(stripes);
		
		Spectra3dPanel panel=new Spectra3dPanel(stripes, ions, params.getFragmentTolerance());
		
		ChartLauncher.openChart(panel.getChart());
	}
}
