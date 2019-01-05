package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import org.jzy3d.chart.ChartLauncher;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import gnu.trove.list.array.TDoubleArrayList;

public class Spectra3dPanelTest {
	private static final int RT_MARGIN=0;
	private static class PeakIntensityComparator implements Comparator<Peak> {
		@Override
		public int compare(Peak o1, Peak o2) {
			if (o1==null&&o2==null) return 0;
			if (o1==null) return -1;
			if (o2==null) return 1;
			return Float.compare(o1.intensity, o2.intensity);
		}
	}

	public static void main(String[] args) throws Exception {
		SearchParameters params=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(16.7), new MassTolerance(16.7), DigestionEnzyme.getEnzyme("trypsin"), DataAcquisitionType.OVERLAPPING_DIA, false, true, false);

		File file=new File("/Users/searleb/Documents/school/perspective/rawfiles/2017dec27_prm_6b_rep1_180103142426.mzML");
		StripeFileInterface raw=StripeFileGenerator.getFile(file, params, true);
		
		float rtInSecStart=71.48f*60;
		float rtInSecStop=72.56f*60;
		
		ArrayList<PrecursorScan> stripes=raw.getPrecursors(rtInSecStart, rtInSecStop);
		Collections.sort(stripes);
		
		ArrayList<Peak> peaks=new ArrayList<>();
		for (PrecursorScan precursorScan : stripes) {
			double[] masses=precursorScan.getMassArray();
			float[] intensities=precursorScan.getIntensityArray();
			for (int i=0; i<intensities.length; i++) {
				if (intensities[i]>10000000) {
					peaks.add(new Peak(masses[i], intensities[i]));
				}
			}
		}
		Collections.sort(peaks, new PeakIntensityComparator());
		
		MassTolerance tolerance=params.getPrecursorTolerance();
		TDoubleArrayList masses=new TDoubleArrayList();
		for (int i=peaks.size()-1; i>=0; i--) {
			double mass=peaks.get(i).mass;
			if (masses.size()==0) {
				masses.add(mass);
			} else {
				int index=masses.indexOf(mass);
				boolean looking=true;
				if (looking&&index>0&&index<=masses.size()) {
					if (tolerance.equals(mass, masses.get(index-1))) looking=false;
				}
				if (looking&&index>=0&&index<masses.size()) {
					if (tolerance.equals(mass, masses.get(index))) looking=false;
				}
				if (looking&&index>=-1&&index<(masses.size()-1)) {
					if (tolerance.equals(mass, masses.get(index+1))) looking=false;
				}
				if (looking) {
					masses.add(mass);
				}
			}
		}
		ArrayList<FragmentIon> ions=new ArrayList<>();
		int index=1;
		for (double mass : masses.toArray()) {
			ions.add(new FragmentIon(mass, (byte)(index), IonType.y));
			index++;
		}
		System.out.println(index+" total peaks");
		
		
		Spectra3dPanel panel=new Spectra3dPanel(stripes, ions.toArray(new FragmentIon[ions.size()]), params.getFragmentTolerance());
		
		ChartLauncher.openChart(panel.getChart());
	}

	public static void main2(String[] args) throws Exception {
		SearchParameters params=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), DataAcquisitionType.OVERLAPPING_DIA, false, true, false);

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
		
		
		ArrayList<FragmentScan> stripes=raw.getStripes(targetMz, rtInSecStart-RT_MARGIN, rtInSecStop+RT_MARGIN, false);
		Collections.sort(stripes);
		
		Spectra3dPanel panel=new Spectra3dPanel(stripes, ions, params.getFragmentTolerance());
		
		ChartLauncher.openChart(panel.getChart());
	}
}
