package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.awt.Dimension;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ScoredPSM;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.ChromatogramExtractor;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Ion;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;

public class ThesaurusTest {
	public static void main(String[] args) throws Exception {
		LibraryFile.OPEN_IN_PLACE=true;
		File libraryFile=new File("/Users/searleb/Documents/phospho_localization/rosenberger/dda_psgs_consensus_spectral_library/psgs_standard_consensus_filtered.dlib");
		File diaFile=new File("/Users/searleb/Documents/phospho_localization/rosenberger/chludwig_K150309_013_SW_0-Pmix_George_dilution_0_in_human.dia");
		LibraryFile library=new LibraryFile();
		library.openFile(libraryFile);
		
		HashMap<String, String> defaults=SearchParameterParser.getDefaultParameters();
		defaults.put("-localizationModification", "Phosphorylation");
		defaults.put("-ptol", "16.67");
		defaults.put("-ftol", "16.67");
		defaults.put("-lftol", "16.67");
		//defaults.put("-frag", "yonly");
		defaults.put("-scoringBreadthType", "uncal20");
		//defaults.put("-scoringBreadthType", "window");
		
		SearchParameters parameters=SearchParameterParser.parseParameters(defaults);
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
		
		UnitBackgroundFrequencyCalculator unitbackground=new UnitBackgroundFrequencyCalculator(0.01f);
		BackgroundFrequencyInterface background=unitbackground;
		//background=BackgroundFrequencyCalculator.generateBackground(stripefile);
		
		PhosphoLocalizer localizer=new PhosphoLocalizer(stripefile, PeptideModification.phosphorylation, background, parameters);

		ArrayList<LibraryEntry> entries=new ArrayList<>();
		LibraryEntry libentry;
		if (true) {
			libentry=library.getEntries("RLS[+79.966331]ISSLNPSNALSR[+10.008269]", (byte)3, false).get(0);
			libentry=libentry.updateRetentionTime(3927.18115234375f);
			entries.add(libentry);
			libentry=library.getEntries("RLSIS[+79.966331]SLNPSNALSR[+10.008269]", (byte)3, false).get(0);
			libentry=libentry.updateRetentionTime(3927.18115234375f);
			entries.add(libentry);
			libentry=library.getEntries("RLSISS[+79.966331]LNPSNALSR[+10.008269]", (byte)3, false).get(0);
			libentry=libentry.updateRetentionTime(3927.18115234375f);
			entries.add(libentry);
		} else {
			System.exit(1);
		}
		
		ArrayList<String> forms=new ArrayList<>();
		HashMap<String, FragmentationModel> models=new HashMap<>();
		for (LibraryEntry entry : entries) {
			models.put(entry.getPeptideModSeq(), PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants()));
			forms.add(entry.getPeptideModSeq());
		}
		
		double precursorMz=entries.get(0).getPrecursorMZ();
		ArrayList<FragmentScan> stripes=stripefile.getStripes(precursorMz, 0, Float.MAX_VALUE, false);
		ArrayList<Spectrum> spectra=FragmentScan.downcastStripeToSpectrum(stripes);
		for (LibraryEntry entry : entries) {
			AnnotatedLibraryEntry anEntry=new AnnotatedLibraryEntry(entry, parameters);

			LibraryEntry bestAlt=null;
			FragmentIon[] localizingIons=null;
			for (LibraryEntry altEntry : entries) {
				if (entry!=altEntry) {
					FragmentIon[] ions=ThesaurusOneScoringTask.getUniqueFragmentIons(models.get(entry.getPeptideModSeq()), models.get(altEntry.getPeptideModSeq()), entry.getPrecursorCharge(), parameters);
					if (localizingIons==null||ions.length<localizingIons.length) {
						// fewer ions separate these
						localizingIons=ions;
						bestAlt=altEntry;
					}
				}
			}
			
			FragmentIon[] ionAnnotations=anEntry.getIonAnnotations();
			ArrayList<FragmentIon> namedIons=new ArrayList<>();
			for (FragmentIon ion : ionAnnotations) {
				if (ion!=null) namedIons.add(ion);
			}
			HashMap<FragmentIon, XYTrace> chromatograms=ChromatogramExtractor.extractFragmentChromatograms(parameters.getFragmentTolerance(), namedIons.toArray(new FragmentIon[namedIons.size()]), spectra, libentry.getRetentionTime(), GraphType.dashedline, true, false);
			chromatograms.putAll(ChromatogramExtractor.extractFragmentChromatograms(parameters.getFragmentTolerance(), localizingIons, spectra, null, GraphType.boldline, true, false));
			Charter.launchChart(entry.getPeptideModSeq()+" Retention Time (min)", "Intensity ("+localizingIons.length+")", false, new Dimension(500, 300), chromatograms.values().toArray(new XYTraceInterface[chromatograms.size()]));
			
			System.out.println(entry.getPeptideModSeq()+" ("+localizingIons.length+") vs "+bestAlt.getPeptideModSeq());
			for (Ion fragmentIon : localizingIons) {
				System.out.println("\t"+fragmentIon+" ("+fragmentIon.getMass()+")");
			}
		}

		ArrayList<Spectrum> limitedPrecursors=new ArrayList<>();
		for (PrecursorScan stripe : stripefile.getPrecursors(libentry.getRetentionTime()-600, libentry.getRetentionTime()+600)) {
			limitedPrecursors.add(stripe);
		}
		Charter.launchChart("Retention Time (min)", "Intensity", false, new Dimension(500, 300), ChromatogramExtractor.extractPrecursorChromatograms(parameters.getPrecursorTolerance(), libentry.getPrecursorMZ(), libentry.getPrecursorCharge(), limitedPrecursors, true, false));
		
		
		BlockingQueue<AbstractScoringResult> resultsQueue=new LinkedBlockingQueue<AbstractScoringResult>();
		BlockingQueue<ModificationLocalizationData> localizationQueue=new LinkedBlockingQueue<ModificationLocalizationData>();

		PrecursorScanMap precursors=new PrecursorScanMap(stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));
		Range range=null;
		for (Range thisRange : stripefile.getRanges().keySet()) {
			if (thisRange.contains(precursorMz)) {
				range=thisRange;
				break;
			}
		}
		float dutyCycle=stripefile.getRanges().get(range).getAverageDutyCycle();
		EncyclopediaOneScorer scorer=new EncyclopediaOneScorer(parameters, unitbackground);
		ThesaurusOneScoringTask task=new ThesaurusOneScoringTask(scorer, entries, stripes, dutyCycle, precursors, localizer, resultsQueue, localizationQueue, parameters);
		//CASiLOneScoringTask task=new CASiLOneScoringTask(scorer, entries, stripes, dutyCycle, precursors, localizer, resultsQueue, localizationQueue, parameters);
		task.call();

		int index=0;
		while (!resultsQueue.isEmpty()) {
			if (!resultsQueue.isEmpty()) {
				AbstractScoringResult result=resultsQueue.take();
				ScoredPSM pair=result.getScoredMSMS();
				index++;
				//for (Pair<ScoredObject<FragmentScan>, float[]> pair : data) {
					System.out.println(index+") "+result.getEntry().getPeptideModSeq()+"\t"+pair.getPrimaryScore()+"\t("+(pair.getMSMS().getScanStartTime()/60f)+" minutes)");
				//}
			} else {
				Thread.sleep(10);
			}
		}
		
		while(!localizationQueue.isEmpty()) {
			if (!localizationQueue.isEmpty()) {
				ModificationLocalizationData data=localizationQueue.take();
				System.out.println(data.getLocalizationPeptideModSeq().getPeptideAnnotation()+" ("+data.isSiteSpecific()+","+data.isLocalized()+") --> "+data.getLocalizationScore()+"\t"+data.getLocalizingIntensity()+"\t"+data.getTotalIntensity()+"\t"+FragmentIon.toArchiveString(data.getLocalizingIons()));
				for (Ion ion : data.getLocalizingIons()) {
					//System.out.println("\t"+ion);
				}
			} else {
				Thread.sleep(10);
			}
		}
	}
}
