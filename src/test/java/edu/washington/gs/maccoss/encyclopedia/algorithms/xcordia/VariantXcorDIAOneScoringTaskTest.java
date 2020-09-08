package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
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
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.UnitBackgroundFrequencyCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.SimilarPeptideBinner;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.VariantFastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideDatabase;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.ChromatogramExtractor;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Ion;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class VariantXcorDIAOneScoringTaskTest {
	public static void main(String[] args) throws Exception {

		HashMap<String, String> defaults=PecanParameterParser.getDefaultParameters();
		defaults.put("-variable", "S=79.966331,T=79.966331,Y=79.966331");
		defaults.put("-localizationModification", "Phosphorylation");
		defaults.put("-requireVariableMods", "true");
		PecanSearchParameters parameters=PecanParameterParser.parseParameters(defaults);
		
		System.out.println("Reading raw file...");
		File diaFile=new File("/Users/searleb/Documents/school/xcordia/phospho/20170430_HeLa_phosp_DIA_B_01_170506220515.dia");
		//File diaFile=new File("/Users/searleb/Documents/xcordia_manuscript/demux/20141121_3_4_DIA_1.dia");
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
		
		UnitBackgroundFrequencyCalculator unitBackgroundFrequencyCalculator=new UnitBackgroundFrequencyCalculator(0.01f);
		UnitBackgroundFrequencyCalculator background=unitBackgroundFrequencyCalculator;
		//background=BackgroundFrequencyCalculator.generateBackground(stripefile);
		
		ArrayList<Range> ranges=new ArrayList<>(stripefile.getRanges().keySet());
		Collections.sort(ranges);

		String proteinFasta=">nxp:NX_O60271-1 \\DbUniqueId=NX_O60271-1 \\PName=C-Jun-amino-terminal kinase-interacting protein 4 isoform Iso 1 \\GName=SPAG9 \\NcbiTaxId=9606 \\TaxName=Homo Sapiens \\Length=1321 \\SV=287 \\EV=597 \\PE=1 \\ModResPsi=(1|MOD:00058|N-acetyl-L-methionine)(109|MOD:00046|O-phospho-L-serine)(177|MOD:00047|O-phospho-L-threonine)(183|MOD:00046|O-phospho-L-serine)(185|MOD:00046|O-phospho-L-serine)(190|MOD:00046|O-phospho-L-serine)(191|MOD:00047|O-phospho-L-threonine)(194|MOD:00046|O-phospho-L-serine)(203|MOD:00046|O-phospho-L-serine)(217|MOD:00047|O-phospho-L-threonine)(226|MOD:00047|O-phospho-L-threonine)(229|MOD:00046|O-phospho-L-serine)(238|MOD:00046|O-phospho-L-serine)(249|MOD:00046|O-phospho-L-serine)(251|MOD:00046|O-phospho-L-serine)(265|MOD:00046|O-phospho-L-serine)(268|MOD:00046|O-phospho-L-serine)(272|MOD:00046|O-phospho-L-serine)(275|MOD:00047|O-phospho-L-threonine)(276|MOD:00047|O-phospho-L-threonine)(279|MOD:00046|O-phospho-L-serine)(280|MOD:00047|O-phospho-L-threonine)(283|MOD:00046|O-phospho-L-serine)(287|MOD:00047|O-phospho-L-threonine)(290|MOD:00047|O-phospho-L-threonine)(292|MOD:00047|O-phospho-L-threonine)(311|MOD:00046|O-phospho-L-serine)(329|MOD:00046|O-phospho-L-serine)(330|MOD:00047|O-phospho-L-threonine)(332|MOD:00046|O-phospho-L-serine)(339|MOD:00046|O-phospho-L-serine)(347|MOD:00046|O-phospho-L-serine)(348|MOD:00047|O-phospho-L-threonine)(363|MOD:00046|O-phospho-L-serine)(364|MOD:00046|O-phospho-L-serine)(365|MOD:00047|O-phospho-L-threonine)(367|MOD:00047|O-phospho-L-threonine)(379|MOD:00047|O-phospho-L-threonine)(381|MOD:00046|O-phospho-L-serine)(387|MOD:00046|O-phospho-L-serine)(388|MOD:00046|O-phospho-L-serine)(391|MOD:00046|O-phospho-L-serine)(418|MOD:00047|O-phospho-L-threonine)(493|MOD:00046|O-phospho-L-serine)(497|MOD:00047|O-phospho-L-threonine)(586|MOD:00047|O-phospho-L-threonine)(588|MOD:00046|O-phospho-L-serine)(593|MOD:00046|O-phospho-L-serine)(594|MOD:00046|O-phospho-L-serine)(595|MOD:00047|O-phospho-L-threonine)(597|MOD:00046|O-phospho-L-serine)(705|MOD:00046|O-phospho-L-serine)(728|MOD:00046|O-phospho-L-serine)(730|MOD:00046|O-phospho-L-serine)(732|MOD:00046|O-phospho-L-serine)(733|MOD:00046|O-phospho-L-serine)(813|MOD:00046|O-phospho-L-serine)(815|MOD:00046|O-phospho-L-serine)(932|MOD:00046|O-phospho-L-serine)(1188|MOD:00046|O-phospho-L-serine)(1238|MOD:00046|O-phospho-L-serine)(1244|MOD:00046|O-phospho-L-serine)(1246|MOD:00047|O-phospho-L-threonine)(1249|MOD:00047|O-phospho-L-threonine)(1262|MOD:00046|O-phospho-L-serine)(1264|MOD:00047|O-phospho-L-threonine)(1302|MOD:00046|O-phospho-L-serine) \\VariantSimple=(8|A)(8|L)(10|E)(12|D)(14|R)(26|D)(30|C)(30|F)(41|C)(41|H)(43|N)(45|V)(59|D)(64|M)(66|T)(69|R)(70|K)(77|V)(80|V)(109|F)(114|R)(115|R)(116|N)(118|H)(120|Q)(121|M)(128|A)(135|V)(137|S)(138|C)(143|G)(145|P)(168|R)(171|V)(180|Y)(185|T)(185|N)(190|A)(192|P)(192|G)(193|Y)(194|N)(200|C)(200|H)(201|S)(203|L)(206|F)(211|D)(212|V)(217|I)(218|A)(223|A)(224|R)(226|S)(227|L)(233|I)(238|C)(241|H)(241|C)(243|D)(243|Y)(247|Q)(256|N)(257|V)(266|V)(274|G)(277|L)(279|L)(281|T)(282|T)(282|D)(283|L)(284|V)(284|G)(285|A)(287|A)(287|R)(291|V)(292|A)(295|N)(299|K)(300|E)(301|S)(302|L)(304|I)(305|I)(306|A)(307|V)(309|H)(311|L)(315|E)(316|R)(317|V)(321|A)(322|D)(322|S)(324|K)(337|Q)(341|F)(342|*)(344|V)(346|K)(349|R)(350|D)(353|V)(356|H)(357|F)(358|N)(366|R)(367|I)(370|V)(376|H)(377|H)(382|V)(386|V)(390|V)(391|L)(393|V)(399|K)(408|Q)(410|A)(412|S)(418|R)(418|A)(424|T)(429|M)(430|A)(431|N)(437|T)(439|G)(458|V)(465|N)(465|R)(467|K)(468|K)(472|Q)(477|Q)(477|W)(480|S)(483|V)(484|K)(487|V)(489|G)(491|N)(492|N)(493|N)(493|C)(495|T)(502|W)(507|K)(521|S)(526|R)(530|*)(533|K)(535|N)(544|I)(545|P)(546|K)(549|S)(554|H)(558|Q)(566|M)(568|T)(569|T)(573|A)(574|F)(579|S)(579|D)(580|V)(580|T)(582|M)(584|Y)(585|I)(586|I)(587|L)(589|I)(589|F)(594|N)(594|G)(604|Y)(606|G)(611|G)(616|G)(616|S)(618|S)(619|V)(621|H)(624|R)(630|C)(630|H)(634|T)(641|S)(646|L)(648|C)(648|*)(651|L)(657|I)(658|A)(661|H)(665|N)(671|M)(673|A)(674|C)(681|T)(684|*)(685|V)(692|A)(693|I)(694|S)(695|V)(695|F)(697|V)(697|S)(699|R)(699|N)(706|I)(709|G)(715|V)(717|V)(718|A)(727|Q)(728|G)(730|C)(733|T)(735|G)(748|S)(749|E)(756|T)(760|V)(762|I)(765|Y)(766|L)(768|I)(770|G)(771|I)(772|T)(779|D)(783|E)(785|L)(786|S)(788|G)(788|*)(789|S)(797|N)(803|V)(805|E)(808|T)(818|A)(820|R)(824|Y)(824|R)(830|S)(831|C)(835|I)(837|I)(840|E)(842|N)(846|D)(847|W)(853|M)(855|P)(856|S)(859|S)(861|I)(49|A)(107|K)(291|N)(557|R)(647|A)(771|F)(863|S)(866|L)(867|L)(868|L)(869|N)(871|R)(872|L)(874|I)(879|R)(881|A)(884|Y)(884|S)(892|A)(895|A)(897|E)(899|V)(901|P)(907|N)(913|I)(915|A)(915|S)(916|*)(917|Y)(920|I)(922|F)(933|A)(934|M)(937|L)(940|G)(949|P)(950|A)(956|N)(958|L)(961|A)(963|L)(967|N)(971|I)(971|A)(972|L)(972|T)(977|H)(979|V)(982|H)(987|I)(987|A)(988|T)(993|F)(997|F)(999|F)(1002|L)(1007|L)(1012|V)(1013|M)(1013|A)(1019|G)(1023|P)(1024|V)(1029|L)(1030|A)(1031|E)(1032|R)(1043|I)(1045|Q)(1045|W)(1049|F)(1050|V)(1051|H)(1054|N)(1054|S)(1057|R)(1063|V)(1069|C)(1070|A)(1076|V)(1076|T)(1083|V)(1083|E)(1084|V)(1085|L)(1086|A)(1093|*)(1096|T)(1096|P)(1097|L)(1102|M)(1107|H)(1111|M)(1112|I)(1114|F)(1116|R)(1122|Q)(1131|C)(1135|L)(1145|C)(1152|I)(1152|V)(1156|S)(1157|H)(1158|F)(1160|L)(1161|E)(1162|I)(1167|V)(1169|C)(1175|I)(1179|P)(1180|D)(1181|I)(1184|K)(1185|H)(1186|S)(1191|H)(1191|C)(1191|P)(1193|C)(1194|V)(1202|S)(1206|C)(1206|S)(1214|Y)(1223|W)(1223|Q)(1225|V)(1228|V)(1228|C)(1230|L)(1237|F)(1241|I)(1241|N)(1244|T)(1246|M)(1253|T)(1254|R)(1255|L)(1259|Q)(1259|K)(1261|D)(1264|M)(1265|L)(1265|T)(1266|F)(1271|A)(1273|N)(1279|V)(1280|N)(1282|Q)(1283|I)(1287|S)(1290|L)(2|Q)(153|N)(165|D)(226|I)(271|R)(322|V)(333|S)(380|D)(385|G)(398|N)(401|S)(577|M)(600|A)(621|C)(675|F)(711|G)(727|L)(799|S)(832|L)(922|L)(961|V)(1003|V)(1005|G)(1086|F)(1087|M)(1098|G)(1144|C)(1150|V)(1174|D)(1192|A)(1196|Q)(1208|L)(1224|N)(1274|R)(1282|P)(1291|A)(1292|F)(1294|R)(1297|I)(1300|K)(1301|L)(1305|N)(1308|M)(1317|I)(1320|S) \\VariantComplex=(269|269|QPQ)(291|291|DD)(337|337|)(563|563|)(565|565|)(753|753|) \\Processed=(1|1321|mature protein)\n" + 
				"MELEDGVVYQEEPGGSGAVMSERVSGLAGSIYREFERLIGRYDEEVVKELMPLVVAVLEN\n" + 
				"LDSVFAQDQEHQVELELLRDDNEQLITQYEREKALRKHAEEKFIEFEDSQEQEKKDLQTR\n" + 
				"VESLESQTRQLELKAKNYADQISRLEEREAELKKEYNALHQRHTEMIHNYMEHLERTKLH\n" + 
				"QLSGSDQLESTAHSRIRKERPISLGIFPLPAGDGLLTPDAQKGGETPGSEQWKFQELSQP\n" + 
				"RSHTSLKVSNSPEPQKAVEQEDELSDVSQGGSKATTPASTANSDVATIPTDTPLKEENEG\n" + 
				"FVKVTDAPNKSEISKHIEVQVAQETRNVSTGSAENEEKSEVQAIIESTPELDMDKDLSGY\n" + 
				"KGSSTPTKGIENKAFDRNTESLFEELSSAGSGLIGDVDEGADLLGMGREVENLILENTQL\n" + 
				"LETKNALNIVKNDLIAKVDELTCEKDVLQGELEAVKQAKLKLEEKNRELEEELRKARAEA\n" + 
				"EDARQKAKDDDDSDIPTAQRKRFTRVEMARVLMERNQYKERLMELQEAVRWTEMIRASRE\n" + 
				"NPAMQEKKRSSIWQFFSRLFSSSSNTTKKPEPPVNLKYNAPTSHVTPSVKKRSSTLSQLP\n" + 
				"GDKSKAFDFLSEETEASLASRREQKREQYRQVKAHVQKEDGRVQAFGWSLPQKYKQVTNG\n" + 
				"QGENKMKNLPVPVYLRPLDEKDTSMKLWCAVGVNLSGGKTRDGGSVVGASVFYKDVAGLD\n" + 
				"TEGSKQRSASQSSLDKLDQELKEQQKELKNQEELSSLVWICTSTHSATKVLIIDAVQPGN\n" + 
				"ILDSFTVCNSHVLCIASVPGARETDYPAGEDLSESGQVDKASLCGSMTSNSSAETDSLLG\n" + 
				"GITVVGCSAEGVTGAATSPSTNGASPVMDKPPEMEAENSEVDENVPTAEEATEATEGNAG\n" + 
				"SAEDTVDISQTGVYTEHVFTDPLGVQIPEDLSPVYQSSNDSDAYKDQISVLPNEQDLVRE\n" + 
				"EAQKMSSLLPTMWLGAQNGCLYVHSSVAQWRKCLHSIKLKDSILSIVHVKGIVLVALADG\n" + 
				"TLAIFHRGVDGQWDLSNYHLLDLGRPHHSIRCMTVVHDKVWCGYRNKIYVVQPKAMKIEK\n" + 
				"SFDAHPRKESQVRQLAWVGDGVWVSIRLDSTLRLYHAHTYQHLQDVDIEPYVSKMLGTGK\n" + 
				"LGFSFVRITALMVSCNRLWVGTGNGVIISIPLTETNKTSGVPGNRPGSVIRVYGDENSDK\n" + 
				"VTPGTFIPYCSMAHAQLCFHGHRDAVKFFVAVPGQVISPQSSSSGTDLTGDKAGPSAQEP\n" + 
				"GSQTPLKSMLVISGGEGYIDFRMGDEGGESELLGEDLPLEPSVTKAERSHLIVWQVMYGN\n" + 
				"E";
		FastaEntryInterface protein=FastaReader.readFasta(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(proteinFasta.getBytes(StandardCharsets.UTF_8)))), "", "", true, parameters).get(0);
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		
		// WITH FIXED AND VARIABLE MODS
		PeptideDatabase targets=new PeptideDatabase();
		/*ArrayList<FastaPeptideEntry> sequences=enzyme.digestProtein(protein, 35, 36, 0, parameters.getAAConstants(), true);
		for (FastaPeptideEntry peptide : sequences) {
			System.out.println(peptide.getSequence());
			targets.add(peptide);
		}*/
		/*targets.add(new FastaPeptideEntry("S[+79.966331]APASPTHPGLMSPR"));
		targets.add(new FastaPeptideEntry("SAPAS[+79.966331]PTHPGLMSPR"));
		targets.add(new FastaPeptideEntry("SAPASPT[+79.966331]HPGLMSPR"));
		targets.add(new FastaPeptideEntry("SAPASPTHPGLMS[+79.966331]PR"));*/

		targets.add(new FastaPeptideEntry("INPY[+79.966331]MSSPCHIEMILTEK"));
		targets.add(new FastaPeptideEntry("INPYMS[+79.966331]SPCHIEMILTEK"));
		targets.add(new FastaPeptideEntry("INPYMSS[+79.966331]PCHIEMILTEK"));
		targets.add(new FastaPeptideEntry("INPYMSSPCHIEMILT[+79.966331]EK"));
		
		System.out.println("Total unique peptides: "+targets.size());

		XCorDIAOneScorer scorer=new XCorDIAOneScorer(parameters, background);
		PrecursorScanMap precursors=new PrecursorScanMap(stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));

		BlockingQueue<PeptideScoringResult> resultsQueue=new LinkedBlockingQueue<PeptideScoringResult>();
		LinkedBlockingQueue<ModificationLocalizationData> localizationQueue=new LinkedBlockingQueue<ModificationLocalizationData>();
		for (Range range : ranges) {
			float dutyCycle=stripefile.getRanges().get(range);
			ArrayList<FragmentScan> stripes=stripefile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, true);
			ArrayList<FragmentScan> xcorStripes=new ArrayList<>();
			for (FragmentScan stripe : stripes) {
				xcorStripes.add(new XCorrStripe(stripe, parameters));
			}
			
			HashSet<FastaPeptideEntry> peptides=XCorDIA.getPeptidesInRange(parameters, targets, range);
			SimilarPeptideBinner binner=new SimilarPeptideBinner();
			ArrayList<ArrayList<FastaPeptideEntry>> bins=binner.binPeptides(peptides);

			for (ArrayList<FastaPeptideEntry> bin : bins) {
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

				LocalizingXcorDIAOneScoringTask task=new LocalizingXcorDIAOneScoringTask(scorer, background, tasks, xcorStripes, range, dutyCycle, precursors, resultsQueue, localizationQueue, parameters);
				//VariantXcorDIAOneScoringTask task=new VariantXcorDIAOneScoringTask(scorer, background, tasks, xcorStripes, range, dutyCycle, precursors, resultsQueue, localizationQueue, parameters);
				task.call();
				
				ArrayList<ModificationLocalizationData> localized=new ArrayList<>();
				float minRT=Float.MAX_VALUE;
				float maxRT=-Float.MAX_VALUE;
				while (localizationQueue.size()>0) {
					ModificationLocalizationData data=localizationQueue.take();
					String peptideModSeq=data.getLocalizationPeptideModSeq().getPeptideModSeq();
					float rtInSeconds=data.getRetentionTimeApexInSeconds();
					Ion[] targetIons=data.getLocalizingIons();
					
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
}
