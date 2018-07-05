package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.XCorDIA;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.BackgroundFrequencyInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.UnitBackgroundFrequencyCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.VariantXcorDIAOneScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorrLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorrStripe;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideDatabase;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.ChromatogramExtractor;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import junit.framework.TestCase;

public class SimilarPeptideBinnerTest extends TestCase {
	public static void main(String[] args) throws Exception {

		HashMap<String, String> defaults=PecanParameterParser.getDefaultParameters();
		defaults.put("-localizationModification", "Phosphorylation");
		defaults.put("-scoringBreadthType", "uncal20");
		PecanSearchParameters parameters=PecanParameterParser.parseParameters(defaults);
		
		System.out.println("Reading raw file...");
		//File diaFile=new File("/Users/searleb/Documents/backup/xcordia_manuscript/xcordia_5p/20141121_27_1_DIA_1.dia");
		File diaFile=new File("/Users/searleb/Documents/xcordia_manuscript/demux/20141121_3_4_DIA_1.dia");
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
		
		UnitBackgroundFrequencyCalculator unitBackgroundFrequencyCalculator=new UnitBackgroundFrequencyCalculator(0.01f);
		UnitBackgroundFrequencyCalculator background=unitBackgroundFrequencyCalculator;
		//background=BackgroundFrequencyCalculator.generateBackground(stripefile);
		
		ArrayList<Range> ranges=new ArrayList<>(stripefile.getRanges().keySet());
		Collections.sort(ranges);

		System.out.println("Reading peff fasta file...");
		File peffFile=new File("/Users/searleb/Documents/xcordia_manuscript/LCM_identified_protein.peff");
		InputStream is=new FileInputStream(peffFile);
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(new BufferedReader(new InputStreamReader(is)), peffFile.getName(), null, true, parameters);
		
		PeptideDatabase targets=new PeptideDatabase();
		for (FastaEntryInterface protein : entries) {
			ArrayList<FastaPeptideEntry> peptides=parameters.getEnzyme().digestProtein(protein, 6, 100, 0, parameters.getAAConstants());
			for (FastaPeptideEntry peptide : peptides) {
				targets.add(peptide);
			}
		}
		
		System.out.println("Total unique peptides: "+targets.size());

		XCorDIAOneScorer scorer=new XCorDIAOneScorer(parameters, background);
		PrecursorScanMap precursors=new PrecursorScanMap(stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));

		BlockingQueue<PeptideScoringResult> resultsQueue=new LinkedBlockingQueue<PeptideScoringResult>();
		LinkedBlockingQueue<ModificationLocalizationData> localizationQueue=new LinkedBlockingQueue<ModificationLocalizationData>();
		for (Range range : ranges) {
			float dutyCycle=stripefile.getRanges().get(range);
			ArrayList<Stripe> stripes=stripefile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, true);
			ArrayList<Stripe> xcorStripes=new ArrayList<>();
			for (Stripe stripe : stripes) {
				xcorStripes.add(new XCorrStripe(stripe, parameters));
			}
			
			HashSet<FastaPeptideEntry> peptides=XCorDIA.getPeptidesInRange(parameters, targets, range);
			SimilarPeptideBinner binner=new SimilarPeptideBinner();
			ArrayList<ArrayList<FastaPeptideEntry>> bins=binner.binPeptides(peptides);

			for (ArrayList<FastaPeptideEntry> bin : bins) {
				boolean keepWorking=true;
				for (FastaPeptideEntry peptide : bin) {
					if (peptide.getAccessions().contains("nxp:NX_P0DJI8-1")) {
						keepWorking=true;
						break;
					}
				}
				if (!keepWorking) continue;

				ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
				HashMap<String, FastaPeptideEntry> entryBySequence=new HashMap<>();
				for (FastaPeptideEntry peptide : bin) {
					for (byte charge=parameters.getMinCharge(); charge<=parameters.getMaxCharge(); charge++) {
						double mz=parameters.getAAConstants().getChargedMass(peptide.getSequence(), charge);
						if (range.contains((float)mz)) {
							XCorrLibraryEntry entry=XCorrLibraryEntry.generateEntry(false, peptide, charge, parameters);
							entryBySequence.put(entry.getPeptideModSeq(), peptide);
							tasks.add(entry);	
						}
					}
				}

				VariantXcorDIAOneScoringTask task=new VariantXcorDIAOneScoringTask(scorer, background, tasks, xcorStripes, range, dutyCycle, precursors, resultsQueue, localizationQueue, parameters);
				task.call();
				
				ArrayList<ModificationLocalizationData> localized=new ArrayList<>();
				float minRT=Float.MAX_VALUE;
				float maxRT=-Float.MAX_VALUE;
				while (localizationQueue.size()>0) {
					ModificationLocalizationData data=localizationQueue.take();
					String peptideModSeq=data.getLocalizationPeptideModSeq().getPeptideModSeq();
					float rtInSeconds=data.getRetentionTimeApexInSeconds();
					FragmentIon[] targetIons=data.getLocalizingIons();
					
					System.out.println(peptideModSeq+"("+targetIons.length+")\trt:"+(rtInSeconds/60.0f)+"\tlocalized:"+data.isLocalized()+"(score:"+data.getLocalizationScore()+")");
					
					//if (data.isLocalized()) {
						localized.add(data);
						if (rtInSeconds>maxRT) maxRT=rtInSeconds;
						if (rtInSeconds<minRT) minRT=rtInSeconds;
					//}
				}
				HashMap<String, ChartPanel> panels=new HashMap<>();
				boolean anyLocalized=false;
				for (ModificationLocalizationData data : localized) {
					if (data.isLocalized()) anyLocalized=true;
					String peptideModSeq=data.getLocalizationPeptideModSeq().getPeptideModSeq();
					for (byte charge=parameters.getMinCharge(); charge<=parameters.getMaxCharge(); charge++) {

						double mz=parameters.getAAConstants().getChargedMass(peptideModSeq, charge);
						if (range.contains((float)mz)) {
							FastaPeptideEntry fastaPeptideEntry = entryBySequence.get(peptideModSeq);
							String variantTag="";
							if (fastaPeptideEntry instanceof VariantFastaPeptideEntry) {
								variantTag=", "+((VariantFastaPeptideEntry) fastaPeptideEntry).getVariant().toString();
							}
							String accessions=General.toString(fastaPeptideEntry.getAccessions());
							float rtInSeconds=data.getRetentionTimeApexInSeconds();
							FragmentIon[] targetIons=data.getLocalizingIons();
							FragmentIon[] allIons=PeptideUtils.getPeptideModel(peptideModSeq, parameters.getAAConstants()).getPrimaryIonObjects(parameters.getFragType(), charge, true);
		
							ArrayList<Spectrum> wideStripeSubset=PhosphoLocalizer.getScanSubsetFromStripes(minRT-60, maxRT+60, stripes);
							HashMap<FragmentIon, XYTrace> uniqueFragmentIons=ChromatogramExtractor.extractFragmentChromatograms(parameters.getFragmentTolerance(), targetIons, wideStripeSubset, (Float)null, GraphType.boldline);
							HashMap<FragmentIon, XYTrace> allFragmentIons=ChromatogramExtractor.extractFragmentChromatograms(parameters.getFragmentTolerance(), allIons, wideStripeSubset, rtInSeconds, GraphType.dashedline);
			
							HashMap<FragmentIon, XYTrace> allFragments=new HashMap<FragmentIon, XYTrace>();
							allFragments.putAll(allFragmentIons);
							allFragments.putAll(uniqueFragmentIons);
							ArrayList<XYTrace> uniqueFragmentsList=new ArrayList<XYTrace>(allFragments.values());
							uniqueFragmentsList.add(new XYTrace(new float[] {rtInSeconds/60f,  rtInSeconds/60f}, new float[] {0.0f, (float)XYTrace.getMaxY(uniqueFragmentsList)}, GraphType.dashedline, "Apex", Color.BLACK, null));
							XYTraceInterface[] fragmentTraces=uniqueFragmentsList.toArray(new XYTrace[uniqueFragmentsList.size()]);
							
							panels.put(peptideModSeq+variantTag, Charter.getChart(accessions+": "+peptideModSeq+" Retention Time (min)", "Intensity", false, fragmentTraces));
						}
					}
				}
				if (anyLocalized&&bin.size()>1&&panels.size()>0) {
					System.out.println("Plotting "+bin.get(0).getSequence());
					Charter.launchCharts(bin.get(0).getSequence(), panels);
				}
			}
		}
		System.out.println("Finished!");
	}
	
	public void testBinner() {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		String annotation=">nxp:NX_P0DJI8-1 \\DbUniqueId=NX_P0DJI8-1 \\PName=Serum amyloid A-1 protein isoform Iso 1 \\GName=SAA1 \\NcbiTaxId=9606 \\TaxName=Homo Sapiens \\Length=122 \\SV=95 \\EV=159 \\PE=1 \\ModResPsi=(101|MOD:00316|N4,N4-dimethyl-L-asparagine) \\VariantSimple=(15|S)(70|A)(75|V)(78|N)(86|L)(90|D) \\Processed=(1|18|signal peptide)(19|94|mature protein)(19|122|mature protein)(20|120|mature protein)(20|121|mature protein)(20|122|mature protein)(21|122|mature protein)(22|119|mature protein)(95|122|maturation peptide)";
		String simpleAnnotation=">nxp:NX_P0DJI8-1";
		String sequence="MKLLTGLVFCSLVLGVSSRSFFSFLGEAFDGARDMWRAYSDMREANYIGSDKYFHARGNYDAAKRGPGGVWAAEAISDARENIQRFFGHGAEDSLADQAANEWGRSGKDPNHFRPAGLPEKY";
		FastaEntry simpleEntry=new FastaEntry("source", simpleAnnotation, sequence);
		ExtendedFastaEntry entry=new ExtendedFastaEntry("source", annotation, sequence, parameters);
		
		ArrayList<FastaPeptideEntry> simplePeptides=parameters.getEnzyme().digestProtein(simpleEntry, 6, 100, 0, parameters.getAAConstants());
		ArrayList<FastaPeptideEntry> peptides=parameters.getEnzyme().digestProtein(entry, 6, 100, 0, parameters.getAAConstants());

		HashSet<FastaPeptideEntry> targets=new HashSet<>();
		for (FastaPeptideEntry peptide : peptides) {
			targets.add(peptide);
		}
		
		SimilarPeptideBinner binner=new SimilarPeptideBinner();
		ArrayList<ArrayList<FastaPeptideEntry>> bins=binner.binPeptides(targets);

		assertFalse(peptides.size()==bins.size());
		assertTrue(simplePeptides.size()==bins.size());
		
		for (ArrayList<FastaPeptideEntry> arrayList : bins) {
			for (FastaPeptideEntry peptide : arrayList) {
				System.out.println(peptide.getSequence());
			}
			System.out.println();
		}
	}
}
