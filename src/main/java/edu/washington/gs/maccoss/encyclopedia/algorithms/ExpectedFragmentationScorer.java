package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import gnu.trove.list.array.TFloatArrayList;

public class ExpectedFragmentationScorer extends AuxillaryPSMScorer {
	private final int startIonIndex=0; //1=start with b2, y2, etc
	

	public ExpectedFragmentationScorer(SearchParameters parameters) {
		super(parameters);
	}

	@Override
	public float[] score(LibraryEntry entry, Stripe spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		FragmentationModel model=new FragmentationModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		MassTolerance tolerance=parameters.getFragmentTolerance();
		double[] masses=spectrum.getMassArray();
		float[] intensities=spectrum.getIntensityArray();
		
		TFloatArrayList ions=new TFloatArrayList();
		FragmentIon[] ys=model.getYIons();
		switch (parameters.getFragType()) {
			case YONLY:
				double[] yIons=FragmentIon.getMasses(ys);
				for (int i=startIonIndex; i<yIons.length; i++) {
					int[] indicies=tolerance.getIndicies(masses, yIons[i]);
					float intensity=0.0f;
					for (int j=0; j<indicies.length; j++) {
						intensity+=intensities[indicies[j]];
					}
					ions.add(intensity);
				}
				if (entry.getPrecursorCharge()>2) {
					yIons=FragmentIon.getMasses(FragmentationModel.getPlus2s(ys));
					for (int i=startIonIndex; i<yIons.length; i++) {
						int[] indicies=tolerance.getIndicies(masses, yIons[i]);
						float intensity=0.0f;
						for (int j=0; j<indicies.length; j++) {
							intensity+=intensities[indicies[j]];
						}
						ions.add(intensity);
					}
				}
				return ions.toArray();
			case CID:
				double[] bIons=FragmentIon.getMasses(model.getBIons());
				for (int i=startIonIndex; i<bIons.length; i++) {
					int[] indicies=tolerance.getIndicies(masses, bIons[i]);
					float intensity=0.0f;
					for (int j=0; j<indicies.length; j++) {
						intensity+=intensities[indicies[j]];
					}
					ions.add(intensity);
				}
				double[] yIonsCID=FragmentIon.getMasses(ys);
				for (int i=startIonIndex; i<yIonsCID.length; i++) {
					int[] indicies=tolerance.getIndicies(masses, yIonsCID[i]);
					float intensity=0.0f;
					for (int j=0; j<indicies.length; j++) {
						intensity+=intensities[indicies[j]];
					}
					ions.add(intensity);
				}
				return ions.toArray();
			case ETD:
				double[] cIons=FragmentIon.getMasses(model.getCIons());
				for (int i=startIonIndex; i<cIons.length; i++) {
					int[] indicies=tolerance.getIndicies(masses, cIons[i]);
					float intensity=0.0f;
					for (int j=0; j<indicies.length; j++) {
						intensity+=intensities[indicies[j]];
					}
					ions.add(intensity);
				}
				double[] zIons=FragmentIon.getMasses(model.getZIons());
				for (int i=startIonIndex; i<zIons.length; i++) {
					int[] indicies=tolerance.getIndicies(masses, zIons[i]);
					float intensity=0.0f;
					for (int j=0; j<indicies.length; j++) {
						intensity+=intensities[indicies[j]];
					}
					ions.add(intensity);
				}
				double[] zp1Ions=FragmentIon.getMasses(model.getZp1Ions());
				for (int i=startIonIndex; i<zp1Ions.length; i++) {
					int[] indicies=tolerance.getIndicies(masses, zp1Ions[i]);
					float intensity=0.0f;
					for (int j=0; j<indicies.length; j++) {
						intensity+=intensities[indicies[j]];
					}
					ions.add(intensity);
				}
				return ions.toArray();
			default:
				throw new EncyclopediaException("Unknown fragmentation type ["+parameters.getFragType()+"]");
		}
	}

	@Override
	public float[] getMissingDataScores(LibraryEntry entry) {
		FragmentationModel model=new FragmentationModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		TFloatArrayList ions=new TFloatArrayList();
		switch (parameters.getFragType()) {
			case YONLY:
				double[] yIons=FragmentIon.getMasses(model.getYIons());
				for (int i=startIonIndex; i<yIons.length; i++) {
					ions.add(0.0f);
				}
				return ions.toArray();
			case CID:
				double[] bIons=FragmentIon.getMasses(model.getBIons());
				for (int i=startIonIndex; i<bIons.length; i++) {
					ions.add(0.0f);
				}
				double[] yIonsCID=FragmentIon.getMasses(model.getYIons());
				for (int i=startIonIndex; i<yIonsCID.length; i++) {
					ions.add(0.0f);
				}
				return ions.toArray();
			case ETD:
				double[] cIons=FragmentIon.getMasses(model.getCIons());
				for (int i=startIonIndex; i<cIons.length; i++) {
					ions.add(0.0f);
				}
				double[] zIons=FragmentIon.getMasses(model.getZIons());
				for (int i=startIonIndex; i<zIons.length; i++) {
					ions.add(0.0f);
				}
				double[] zp1Ions=FragmentIon.getMasses(model.getZp1Ions());
				for (int i=startIonIndex; i<zp1Ions.length; i++) {
					ions.add(0.0f);
				}
				return ions.toArray();
			default:
				throw new EncyclopediaException("Unknown fragmentation type ["+parameters.getFragType()+"]");
		}
	}

	@Override
	public String[] getScoreNames(LibraryEntry entry) {
		FragmentationModel model=new FragmentationModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		ArrayList<String> names=new ArrayList<String>();
		switch (parameters.getFragType()) {
			case YONLY:
				FragmentIon[] yIons=model.getYIons();
				for (int i=startIonIndex; i<yIons.length; i++) {
					names.add("y"+(i+1));
				}
				if (entry.getPrecursorCharge()>2) {
					for (int i=startIonIndex; i<yIons.length; i++) {
						names.add("y"+(i+1)+"+2H");
					}
				}
				return names.toArray(new String[names.size()]);
			case CID:
				FragmentIon[] bIons=model.getBIons();
				for (int i=startIonIndex; i<bIons.length; i++) {
					names.add("b"+(i+1));
				}
				FragmentIon[] yIonsCID=model.getYIons();
				for (int i=startIonIndex; i<yIonsCID.length; i++) {
					names.add("y"+(i+1));
				}
				return names.toArray(new String[names.size()]);
			case ETD:
				FragmentIon[] cIons=model.getCIons();
				for (int i=startIonIndex; i<cIons.length; i++) {
					names.add("c"+(i+1));
				}
				FragmentIon[] zIons=model.getZIons();
				for (int i=startIonIndex; i<zIons.length; i++) {
					names.add("z"+(i+1));
				}
				FragmentIon[] zp1Ions=model.getZp1Ions();
				for (int i=startIonIndex; i<zp1Ions.length; i++) {
					names.add("z"+(i+1)+"+1");
				}
				return names.toArray(new String[names.size()]);
			default:
				throw new EncyclopediaException("Unknown fragmentation type ["+parameters.getFragType()+"]");
		}
	}
}
