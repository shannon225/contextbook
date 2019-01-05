package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.XCorDIA;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.UnitBackgroundFrequencyCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.VariantXcorDIAOneScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorrLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorrStripe;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
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
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import junit.framework.TestCase;

public class SimilarPeptideBinnerTest extends TestCase {
	
	public static void main(String[] args) throws Exception {
		HashMap<String, String> defaults=PecanParameterParser.getDefaultParameters();
		PecanSearchParameters parameters=PecanParameterParser.parseParameters(defaults);

		File f=new File("/Users/searleb/Documents/school/uniprot-9606.fasta");
		//File f=new File("/Users/searleb/Documents/school/xcordia/nextprot2017_testPEFF1.0rc25_a.peff");
		//File f=new File("/Users/searleb/Documents/school/xcordia/LCM_identified_protein.peff");
		InputStream is=new FileInputStream(f);
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(new BufferedReader(new InputStreamReader(is)), f.getName(), null, true, parameters);

		System.out.println("Number of entries: "+entries.size());
		PeptideDatabase targets=new PeptideDatabase();
		for (FastaEntryInterface protein : entries) {
			ArrayList<FastaPeptideEntry> peptides=parameters.getEnzyme().digestProtein(protein, 6, 50, 0, parameters.getAAConstants(), false);
			for (FastaPeptideEntry peptide : peptides) {
				targets.add(peptide);
			}
		}
		System.out.println("Number of peptides: "+targets.size());

		int singletons=0;
		int multiples=0;

		ArrayList<Range> ranges=new ArrayList<>();
		File[] diaFiles=new File[] {new File("/Users/searleb/Documents/school/xcordia/hela/122715_bcs_hela_24mz_400_1000.dia")};
		diaFiles=new File[] {new File("/Users/searleb/Documents/school/xcordia/hela_prm/2017dec27_pool_library1.dia"),
				new File("/Users/searleb/Documents/school/xcordia/hela_prm/2017dec27_pool_library2.dia"),
				new File("/Users/searleb/Documents/school/xcordia/hela_prm/2017dec27_pool_library3.dia"),
				new File("/Users/searleb/Documents/school/xcordia/hela_prm/2017dec27_pool_library4.dia"),
				new File("/Users/searleb/Documents/school/xcordia/hela_prm/2017dec27_pool_library5.dia"),
				new File("/Users/searleb/Documents/school/xcordia/hela_prm/2017dec27_pool_library6.dia")};
		diaFiles=new File[] {new File("/Users/searleb/Documents/school/xcordia/hela_prm/2017dec27_normal_dia_6b_rep1.dia")};
		for (File file : diaFiles) {
			StripeFileInterface stripefile=StripeFileGenerator.getFile(file, parameters);
			ranges.addAll(stripefile.getRanges().keySet());
		}
		Collections.sort(ranges);
		System.out.println("Number of ranges: "+ranges.size());
		
		for (Range range : ranges) {
			System.out.print(range.toString());
			HashSet<FastaPeptideEntry> peptides=XCorDIA.getPeptidesInRange(parameters, targets, range);
			System.out.print("\tinrange:"+peptides.size());
			SimilarPeptideBinner binner=new SimilarPeptideBinner();
			ArrayList<ArrayList<FastaPeptideEntry>> bins=binner.binPeptides(peptides);
			System.out.print("\ttotal:"+targets.size());

			int localsingletons=0;
			int localmultiples=0;
			for (ArrayList<FastaPeptideEntry> list : bins) {
				if (list.size()>1) {
					localmultiples+=list.size();
				} else if (list.size()==1) {
					localsingletons++;
				}
			}
			System.out.println("\tsingle:"+localsingletons+"\tmultiple:"+localmultiples);
			singletons+=localsingletons;
			multiples+=localmultiples;
		}
		System.out.println("FINAL: single:"+singletons+"\tmultiple:"+multiples+"\t"+targets.size());
	}
	public static void main2(String[] args) throws Exception {

		HashMap<String, String> defaults=PecanParameterParser.getDefaultParameters();
		defaults.put("-localizationModification", "Phosphorylation");
		PecanSearchParameters parameters=PecanParameterParser.parseParameters(defaults);
		
		System.out.println("Reading raw file...");
		File diaFile=new File("/Users/searleb/Documents/backup/xcordia_manuscript/xcordia_5p/20141121_27_1_DIA_1.dia");
		//File diaFile=new File("/Users/searleb/Documents/xcordia_manuscript/demux/20141121_3_4_DIA_1.dia");
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
		
		UnitBackgroundFrequencyCalculator unitBackgroundFrequencyCalculator=new UnitBackgroundFrequencyCalculator(0.01f);
		UnitBackgroundFrequencyCalculator background=unitBackgroundFrequencyCalculator;
		//background=BackgroundFrequencyCalculator.generateBackground(stripefile);
		
		ArrayList<Range> ranges=new ArrayList<>(stripefile.getRanges().keySet());
		Collections.sort(ranges);

		System.out.println("Reading peff fasta file...");
		File peffFile=new File("/Users/searleb/Documents/school/xcordia/LCM_identified_protein.peff");
		//File peffFile=new File("/Users/searleb/Documents/xcordia_manuscript/LCM_identified_protein.peff");
		InputStream is=new FileInputStream(peffFile);
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(new BufferedReader(new InputStreamReader(is)), peffFile.getName(), null, true, parameters);
		
		PeptideDatabase targets=new PeptideDatabase();
		for (FastaEntryInterface protein : entries) {
			ArrayList<FastaPeptideEntry> peptides=parameters.getEnzyme().digestProtein(protein, 6, 100, 0, parameters.getAAConstants(), false);
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
		
		ArrayList<FastaPeptideEntry> simplePeptides=parameters.getEnzyme().digestProtein(simpleEntry, 6, 100, 0, parameters.getAAConstants(), false);
		ArrayList<FastaPeptideEntry> peptides=parameters.getEnzyme().digestProtein(entry, 6, 100, 0, parameters.getAAConstants(), false);

		HashSet<FastaPeptideEntry> targets=new HashSet<>();
		for (FastaPeptideEntry peptide : peptides) {
			targets.add(peptide);
		}
		
		SimilarPeptideBinner binner=new SimilarPeptideBinner();
		ArrayList<ArrayList<FastaPeptideEntry>> bins=binner.binPeptides(targets);

		assertFalse(peptides.size()==bins.size());
		assertTrue(simplePeptides.size()==bins.size());
	}
	
	public void testModBinner() {
		String fakeTTR=">nxp:NX_Q96B36-1 \\DbUniqueId=NX_Q96B36-1 \\PName=Proline-rich AKT1 substrate 1 isoform Iso 1 \\GName=AKT1S1 \\NcbiTaxId=9606 \\TaxName=Homo Sapiens \\Length=256 \\SV=135 \\EV=428 \\PE=1 "
				+ "\\ModResPsi=(3|MOD:00046|O-phospho-L-serine)(51|MOD:00078|omega-N-methyl-L-arginine)(73|MOD:00047|O-phospho-L-threonine)(90|MOD:00047|O-phospho-L-threonine)(97|MOD:00047|O-phospho-L-threonine)(116|MOD:00046|O-phospho-L-serine)(187|MOD:00046|O-phospho-L-serine)(198|MOD:00047|O-phospho-L-threonine)(247|MOD:00046|O-phospho-L-serine)(88|MOD:00046|O-phospho-L-serine)(92|MOD:00046|O-phospho-L-serine)(183|MOD:00046|O-phospho-L-serine)(202|MOD:00046|O-phospho-L-serine)(203|MOD:00046|O-phospho-L-serine)(211|MOD:00046|O-phospho-L-serine)(212|MOD:00046|O-phospho-L-serine)(246|MOD:00047|O-phospho-L-threonine) "
				//+ "\\VariantSimple=(1|V)(3|L)(5|H)(12|T)(13|M)(19|C)(21|W)(21|Q)(23|Q)(26|M)(27|D)(33|T)(40|H)(41|L)(42|S)(43|S)(46|H)(46|C)(47|P)(51|Q)(51|L)(55|V)(59|H)(60|C)(60|H)(61|Y)(63|R)(68|G)(76|Q)(76|W)(78|A)(79|V)(80|S)(82|Q)(86|L)(87|S)(88|C)(88|R)(93|T)(95|W)(95|Q)(95|P)(97|A)(99|V)(103|D)(106|H)(107|K)(114|G)(115|I)(122|T)(122|V)(129|I)(130|A)(135|T)(136|I)(141|H)(141|S)(142|R)(142|H)(142|S)(145|K)(147|E)(149|K)(151|I)(161|S)(162|T)(163|D)(166|I)(168|L)(170|H)(172|V)(174|T)(182|E)(185|S)(186|L)(188|L)(190|F)(190|I)(192|S)(197|T)(198|I)(201|W)(201|Q)(204|N)(207|D)(209|R)(209|L)(214|Y)(215|P)(223|H)(233|A)(233|S)(234|H)(239|V)(242|L)(245|K)(117|R)(212|L)(249|L)(252|Q)(7|K)(51|*)(53|S)(58|V)(64|N)(104|K)(158|*)(163|A)(199|*)(200|V) \\VariantComplex=(39|39|PP)(206|206|)(114|115|) \\Processed=(1|256|mature protein)"
				+"\n"+"MASGRPEELWEAVVGAAERFRARTGTELVLLTAAPPPPPRPGPCAYAAHGRGALAEAARR" + 
				"CLHDIALAHRAATAARPPAPPPAPQPPSPTPSPPRPTLAREDNEEDEDEPTETETSGEQL" + 
				"GISDNGGLFVMDEDATLQDLPPFCESDPESTDDGSLSEETPAGPPTCSVPPASALPTQQY" + 
				"AKSLPVSVPVWGFKEKRTEARSSDEENGPPSSPDLDRIAASMRALVLREAEDTQVFGDLP" + 
				"RPRLNTSDFQKLKRKY";
		HashMap<String, String> paramMap=SearchParameterParser.getDefaultParameters();
		paramMap.put("-variable", "S=79.966331,T=79.966331,Y=79.966331");
		SearchParameters parameters=PecanParameterParser.parseParameters(paramMap);
		FastaEntryInterface entry=FastaReader.readFasta(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(fakeTTR.getBytes(StandardCharsets.UTF_8)))), "", "", true, parameters).get(0);
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		
		// WITH FIXED AND VARIABLE MODS
		ArrayList<FastaPeptideEntry> targets=enzyme.digestProtein(entry, 4, 40, 0, parameters.getAAConstants(), false);
		
		SimilarPeptideBinner binner=new SimilarPeptideBinner();
		ArrayList<ArrayList<FastaPeptideEntry>> bins=binner.binPeptides(targets);
		
		assertEquals(12, bins.size());
		
		int count=0;
		for (ArrayList<FastaPeptideEntry> arrayList : bins) {
			count+=arrayList.size();
			for (FastaPeptideEntry peptide : arrayList) {
				System.out.println(peptide.getSequence());
			}
			System.out.println();
		}
		assertEquals(22, count);
	}
}
