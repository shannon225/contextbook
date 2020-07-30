package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class LibraryEntryModifier {
	public static LibraryEntry modifyModelAtEverySite(LibraryEntry entry, TCharDoubleHashMap fixedMods, boolean changePTMs, SearchParameters parameters) {
		FragmentationModel model=modifyModelAtEverySite(entry.getPeptideModSeq(), fixedMods, changePTMs, parameters);
		return entry.getEntryFromNewSequence(model.getPeptideModSeq(), entry.getAccessions(), false, parameters).y;
	}
	public static ArrayList<LibraryEntry> modifyModelAtEachSite(LibraryEntry entry, TCharDoubleHashMap fixedMods, boolean changePTMs, SearchParameters parameters) {
		ArrayList<LibraryEntry> entries=new ArrayList<>();
		ArrayList<FragmentationModel> models=modifyModelAtEachSite(entry.getPeptideModSeq(), fixedMods, changePTMs, parameters);
		for (FragmentationModel model : models) {
			entries.add(entry.getEntryFromNewSequence(model.getPeptideModSeq(), entry.getAccessions(), false, parameters).y);
		}
		return entries;
	}

	/**
	 * 
	 * @param peptideModSeq
	 * @param fixedMods
	 * @param changePTMs true if remove old PTMs at AA and then add new mass, false if combine masses (e.g. SIL)
	 * @param parameters
	 * @return
	 */
	public static FragmentationModel modifyModelAtEverySite(String peptideModSeq, TCharDoubleHashMap fixedMods, boolean changePTMs, SearchParameters parameters) {
		FragmentationModel model=PeptideUtils.getPeptideModel(peptideModSeq, parameters.getAAConstants());
		double[] neutralLosses=model.getNeutralLosses().clone();
		double[] modificationMasses=model.getModificationMasses().clone();
		double[] masses=model.getMasses().clone();
		String[] aas=model.getAas().clone();
		
		for (int i = 0; i < aas.length; i++) {
			char aa=aas[i].charAt(0);
			double ptmMass=fixedMods.get(aa);
			
			if (ptmMass!=0.0) {
				if (changePTMs) {
					aas[i]=aa+"["+ptmMass+"]";
					masses[i]=masses[i]-modificationMasses[i]+ptmMass;
					modificationMasses[i]=ptmMass;
					neutralLosses[i]=0.0; // don't adjust or keep neutral loss intensities since they most likely won't equate
				} else {
					masses[i]=masses[i]+ptmMass;
					modificationMasses[i]=modificationMasses[i]+ptmMass;
					aas[i]=aa+"["+modificationMasses[i]+"]";
					// keep original neutral loss, since the new fixed mass change affects the AA, not the original ptm 
				}
			}
		}
		return new FragmentationModel(masses, modificationMasses, neutralLosses, aas);
	}

	/**
	 * 
	 * @param peptideModSeq
	 * @param fixedMods
	 * @param changePTMs true if remove old PTMs at AA and then add new mass, false if combine masses (e.g. SIL)
	 * @param parameters
	 * @return
	 */
	public static ArrayList<FragmentationModel> modifyModelAtEachSite(String peptideModSeq, TCharDoubleHashMap fixedMods, boolean changePTMs, SearchParameters parameters) {
		FragmentationModel model=PeptideUtils.getPeptideModel(peptideModSeq, parameters.getAAConstants());
		ArrayList<FragmentationModel> modifiedModels=new ArrayList<>();

		String[] originalAAs=model.getAas();
		for (int i = 0; i < originalAAs.length; i++) {
			char aa=originalAAs[i].charAt(0);
			double ptmMass=fixedMods.get(aa);

			if (ptmMass!=0.0) {
				double[] neutralLosses=model.getNeutralLosses().clone();
				double[] modificationMasses=model.getModificationMasses().clone();
				double[] masses=model.getMasses().clone();
				String[] aas=model.getAas().clone();
				if (changePTMs) {
					aas[i]=aa+"["+ptmMass+"]";
					masses[i]=masses[i]-modificationMasses[i]+ptmMass;
					modificationMasses[i]=ptmMass;
					neutralLosses[i]=0.0; // don't adjust or keep neutral loss intensities since they most likely won't equate
				} else {
					masses[i]=masses[i]+ptmMass;
					modificationMasses[i]=modificationMasses[i]+ptmMass;
					aas[i]=aa+"["+modificationMasses[i]+"]";
					// keep original neutral loss, since the new fixed mass change affects the AA, not the original ptm 
				}
				
				FragmentationModel newModel=new FragmentationModel(masses, modificationMasses, neutralLosses, aas);
				modifiedModels.add(newModel);
			}
		}
		return modifiedModels;
	}
}
