package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.io.File;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.LibraryUtilities;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class LibraryEntryModifier {
	public static void main(String[] args) throws Exception {
		//File inputFile=new File("/Users/searleb/Downloads/older_downloads/subset_yeast_dda_library.dlib");
		//File outputFile=new File("/Users/searleb/Downloads/older_downloads/cys_subset_yeast_dda_library.dlib");

		TCharDoubleHashMap ptms=new TCharDoubleHashMap();
		//ptms.put('K',8.014199);
		//ptms.put('R',10.008269);
		
		//ptms.put('C',4.02);
		
		//ptms.put('K',4.025107);
		//ptms.put('R',6.020129);
		
		//LibraryInterface library=BlibToLibraryConverter.getFile(inputFile);
		//LibraryUtilities.modifyLibrary(outputFile, ptms, false, library);
		
		// VARIABLE TESTS
//		File inputFile=new File("/Users/searleb/Downloads/pan_human_library.dlib");
//		File outputFile=new File("/Users/searleb/Downloads/pan_human_library_p12.dlib");
//		ptms.put('n', 12.0);
		
		File inputFile=new File("/Users/searleb/Documents/iarpa/IARPA_bone/10p/temp_bone_refs_vars.dlib");
		File outputFile=new File("/Users/searleb/Documents/iarpa/IARPA_bone/10p/temp_bone_refs_vars_hydroxyproline.dlib");
		ptms.put('P', 15.9949);
		
		LibraryInterface library=BlibToLibraryConverter.getFile(inputFile);
		LibraryUtilities.modifyLibrary(outputFile, ptms, false, library);
	}
	
	public static LibraryEntry modifyModelAtEverySite(LibraryEntry entry, TCharDoubleHashMap fixedMods, boolean changePTMs, SearchParameters parameters) {
		FragmentationModel model=modifyModelAtEverySite(entry.getPeptideModSeq(), fixedMods, changePTMs, parameters);
		return entry.getEntryFromNewSequence(model.getPeptideModSeq(), false, false, false, parameters).y;
	}
	public static ArrayList<LibraryEntry> modifyModelAtEachSite(LibraryEntry entry, TCharDoubleHashMap fixedMods, boolean changePTMs, SearchParameters parameters) {
		ArrayList<LibraryEntry> entries=new ArrayList<>();
		ArrayList<FragmentationModel> models=modifyModelAtEachSite(entry.getPeptideModSeq(), fixedMods, changePTMs, parameters);
		for (FragmentationModel model : models) {
			entries.add(entry.getEntryFromNewSequence(model.getPeptideModSeq(), false, false, false, parameters).y);
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
			checkAndAddPTM(changePTMs, model, modifiedModels, originalAAs, i, ptmMass);
		}

		checkAndAddPTM(changePTMs, model, modifiedModels, originalAAs, 0, fixedMods.get(AminoAcidConstants.N_TERM));
		checkAndAddPTM(changePTMs, model, modifiedModels, originalAAs, originalAAs.length-1, fixedMods.get(AminoAcidConstants.C_TERM));
		
		return modifiedModels;
	}

	private static void checkAndAddPTM(boolean changePTMs, FragmentationModel model,
			ArrayList<FragmentationModel> modifiedModels, String[] originalAAs, int i, double ptmMass) {
		if (ptmMass!=0.0) {
			char aaChar=originalAAs[i].charAt(0);
			double[] neutralLosses=model.getNeutralLosses().clone();
			double[] modificationMasses=model.getModificationMasses().clone();
			double[] masses=model.getMasses().clone();
			String[] aas=model.getAas().clone();
			if (changePTMs) {
				aas[i]=aaChar+"["+ptmMass+"]";
				masses[i]=masses[i]-modificationMasses[i]+ptmMass;
				modificationMasses[i]=ptmMass;
				neutralLosses[i]=0.0; // don't adjust or keep neutral loss intensities since they most likely won't equate
			} else {
				masses[i]=masses[i]+ptmMass;
				modificationMasses[i]=modificationMasses[i]+ptmMass;
				aas[i]=aaChar+"["+modificationMasses[i]+"]";
				// keep original neutral loss, since the new fixed mass change affects the AA, not the original ptm 
			}
			
			FragmentationModel newModel=new FragmentationModel(masses, modificationMasses, neutralLosses, aas);
			modifiedModels.add(newModel);
		}
	}
}
