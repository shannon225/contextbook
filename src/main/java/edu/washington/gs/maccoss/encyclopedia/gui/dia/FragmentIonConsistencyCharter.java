package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.awt.Color;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Optional;

import org.jfree.chart.ChartPanel;
import org.jfree.data.category.DefaultCategoryDataset;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.ChromatogramExtractor;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Correlation;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class FragmentIonConsistencyCharter {
	public static final DecimalFormat df = new DecimalFormat( "#,###,###,##0.00" );

	public static ChartPanel getBarChart(PeptidePrecursor peptide, Optional<Spectrum> libraryEntry, Spectrum[] spectra, String[] sampleName, SearchParameters parameters) {
		assert(spectra.length>0);
		assert(spectra.length==sampleName.length);
		
		DefaultCategoryDataset result = new DefaultCategoryDataset();
		
		FragmentationModel model=PeptideUtils.getPeptideModel(peptide.getPeptideModSeq(), parameters.getAAConstants());
		FragmentIon[] primaryIons = model.getPrimaryIonObjects(parameters.getFragType(), peptide.getPrecursorCharge(), true);
		
		float[] libraryAnnotatedIntensities=new float[primaryIons.length];
		if (libraryEntry.isPresent()) {
			Spectrum entry=libraryEntry.get();
			float[] intensities=entry.getIntensityArray();
			double[] masses=entry.getMassArray();
			
			for (int j = 0; j < primaryIons.length; j++) {
				int[] indicies=parameters.getFragmentTolerance().getIndicies(masses, primaryIons[j].getMass());
				float totalIntensity=0.0f;
				for (int k = 0; k < indicies.length; k++) {
					totalIntensity+=intensities[indicies[k]];
				}
				libraryAnnotatedIntensities[j]=totalIntensity;
			}
		}
		
		HashSet<FragmentIon> foundIons=new HashSet<>();
		
		for (int i = 0; i < spectra.length; i++) {
			if (spectra[i]==null) {
				for (int j = 0; j < primaryIons.length; j++) {
					result.addValue(0.0, primaryIons[j], sampleName[i]);
				}
			} else {
				float[] intensities=spectra[i].getIntensityArray();
				double[] masses=spectra[i].getMassArray();
				
				float[] annotatedIntensities=new float[primaryIons.length];
				for (int j = 0; j < primaryIons.length; j++) {
					int[] indicies=parameters.getFragmentTolerance().getIndicies(masses, primaryIons[j].getMass());
					float totalIntensity=0.0f;
					for (int k = 0; k < indicies.length; k++) {
						totalIntensity+=intensities[indicies[k]];
					}
					annotatedIntensities[j]=totalIntensity;
				}
				
				float totalIntensity=General.sum(annotatedIntensities);

				String name=sampleName[i];
				if (libraryEntry.isPresent()) {
					name=sampleName[i]+" ("+df.format(Correlation.getPearsons(libraryAnnotatedIntensities, annotatedIntensities))+")";
				}
				for (int j = 0; j < annotatedIntensities.length; j++) {
					if (annotatedIntensities[j]>0.0f) {
						foundIons.add(primaryIons[j]);
						result.addValue(annotatedIntensities[j]/totalIntensity, primaryIons[j], name);
					}
				}
			}
		}

		if (libraryEntry.isPresent()) {
			Spectrum entry=libraryEntry.get();
			String name="Library";
			
			float[] intensities=entry.getIntensityArray();
			double[] masses=entry.getMassArray();
			
			FragmentIon[] foundIonArray=foundIons.toArray(new FragmentIon[foundIons.size()]);
			float[] annotatedIntensities=new float[foundIonArray.length];
			for (int j = 0; j < foundIonArray.length; j++) {
				int[] indicies=parameters.getFragmentTolerance().getIndicies(masses, foundIonArray[j].getMass());
				float totalIntensity=0.0f;
				for (int k = 0; k < indicies.length; k++) {
					totalIntensity+=intensities[indicies[k]];
				}
				annotatedIntensities[j]=totalIntensity;
			}
			
			float totalIntensity=General.sum(annotatedIntensities);
			for (int j = 0; j < annotatedIntensities.length; j++) {
				result.addValue(annotatedIntensities[j]/totalIntensity, foundIonArray[j], name);
			}
		}
		
		return Charter.getBarChart(null, "Sample", "Fractional Intensity", result, true);
	}
	
	public static ArrayList<XYTrace> getPrecursorButterfly(byte[] isotopes, float[] top, float[] bottom, Color[] colors) {
		System.out.println(General.toString(top));
		
		float max = General.max(top);
		if (max>0) top=General.divide(top, max);
		max = General.max(bottom);
		if (max>0) bottom=General.divide(bottom, max);
		
		ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
		for (int i = 0; i < isotopes.length; i++) {
			GraphType lineType=isotopes[i]<0?GraphType.bolddashedline:GraphType.boldline;
			
			traces.add(new XYTrace(new double[] {isotopes[i], isotopes[i]}, new double[] {-bottom[i], top[i]}, lineType, getIsotopeShortName(isotopes[i]), colors[i], 3f));
			traces.add(new XYTrace(new double[] {isotopes[i], isotopes[i]}, new double[] {-bottom[i], top[i]}, GraphType.bigpoint, getIsotopeShortName(isotopes[i]), colors[i], 10f));				
		}
		return traces;
	}

	private static String getIsotopeShortName(byte isotope) {
		String name;
		if (isotope>0) {
			name="P+"+isotope;
		} else if (isotope<0) {
			name="P"+isotope;
		} else {
			name="P";
		}
		return name;
	}

	public static LibraryEntry getButterfly(LibraryEntry top, LibraryEntry bottom) {
		double[] masses = new double[top.getMassArray().length + bottom.getMassArray().length];
		System.arraycopy(top.getMassArray(), 0, masses, 0, top.getMassArray().length);
		System.arraycopy(bottom.getMassArray(), 0, masses, top.getMassArray().length, bottom.getMassArray().length);
		float[] intensities = new float[top.getIntensityArray().length + bottom.getIntensityArray().length];
		System.arraycopy(normalize(top.getIntensityArray()), 0, intensities, 0, top.getIntensityArray().length);
		System.arraycopy(General.multiply(normalize(bottom.getIntensityArray()), -1f), 0, intensities,
				top.getIntensityArray().length, bottom.getIntensityArray().length);
		float[] correlations = new float[top.getCorrelationArray().length + bottom.getCorrelationArray().length];
		System.arraycopy(top.getCorrelationArray(), 0, correlations, 0, top.getCorrelationArray().length);
		System.arraycopy(bottom.getCorrelationArray(), 0, correlations, top.getCorrelationArray().length,
				bottom.getCorrelationArray().length);
		boolean[] quantifiedIons = new boolean[top.getQuantifiedIonsArray().length
				+ bottom.getQuantifiedIonsArray().length];
		System.arraycopy(top.getQuantifiedIonsArray(), 0, quantifiedIons, 0, top.getQuantifiedIonsArray().length);
		System.arraycopy(bottom.getQuantifiedIonsArray(), 0, quantifiedIons, top.getQuantifiedIonsArray().length,
				bottom.getQuantifiedIonsArray().length);

		LibraryEntry trace = new LibraryEntry(top.getSource(), top.getAccessions(), top.getSpectrumIndex(),
				top.getPrecursorMZ(), top.getPrecursorCharge(), top.getLegacyPeptideModSeq(), top.getPeptideModSeq(),
				top.getCopies(), top.getRetentionTime(), top.getScore(), masses, intensities, correlations,
				quantifiedIons, top.getIonMobility(), true);
		return trace;
	}

	private static float[] normalize(float[] intensities) {
		return General.divide(intensities, General.max(intensities));
	}
}
