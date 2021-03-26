package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Optional;

import org.jfree.chart.ChartPanel;
import org.jfree.data.category.DefaultCategoryDataset;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
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
}
