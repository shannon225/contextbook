package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;

import com.google.common.util.concurrent.ThreadFactoryBuilder;

import edu.washington.gs.maccoss.encyclopedia.algorithms.BackgroundGenerator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.DotProduct;
import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PecanAuxillaryScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PecanRawScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PecanScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringTask;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PecanLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TDoubleIntHashMap;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.procedure.TDoubleObjectProcedure;
import gnu.trove.set.hash.TDoubleHashSet;

public class StripeExtractor {
	private static final SearchParameters PARAMETERS=new SearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"));

	public static void main(String[] args) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		int cores=Runtime.getRuntime().availableProcessors();

		File f=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/82593_lv_mcx_DIA_5mz_400to525.dia");
		StripeFile stripefile=new StripeFile();
		stripefile.openFile(f);

		// File lf=new
		// File("/Users/searleb/Documents/school/projects/qe_phospho.elib");
		// LibraryFile libraryFile=new LibraryFile();
		// libraryFile.openFile(lf);

		byte[] charges=new byte[] {(byte)2, (byte)3};
		String[] peptides=new String[] {"ADVDAATLAR", "APIQWEER", "ATNLTVSAVR", "AVDSLVPIGR", "AYIDSTDSR", "DGLTDVYNK", "DGPGFYTTR", "DTPGFIVNR", "DVLSNLIPK", "DYPLIPVGK", "EAFLLFDR", "EALISQLTR",
				"FGGGSVELLK", "GGASDALLYR", "GIWHNYDK", "GLNEAAIVNK", "GPLVQGVDSR", "GSGIQWDLR", "GTAVVNGEFK", "IGGIGTVPVGR", "IGYPAPNFK", "ILQEGVDPK", "IRC[+57]DIANVK", "ITQSNAILR", "IVVHAGGVIR",
				"KFVADGIFK", "LGFMSAFVK", "LLEAASVSSK", "LLFEELVR", "LNVLANVIR", "LQGDLVTIR", "LTLSALIDGK", "LTLSALVDGK", "LVNMLDAVR", "MFASFPTTK", "NDAGYSEPR", "NTYYASIAK", "QGVLTLEIR", "QIFLGGVDR",
				"RAEVLDSTK", "RGDFIPGLR", "RLTDADAMK", "RLVVQQAGK", "RWEVAALR", "SFLPLLRR", "SLHTLFGDK", "SVQAAMEKR", "TAYVGENVR", "TDLTAVPASR", "TIPWLENR", "TLEDILFR", "TQLVSNLKK", "TQVQSVIDK",
				"TSGGAGGLGSLR", "VAAENQYGR", "VDFDDIHR", "VGPANPSLQK", "VLSIGDGIAR", "VTLVSAAPEK", "VVDDELATR", "VVFIFGPDK", "WPLYLSTK", "YDHLGDSPK", "YVDMSAKSK", "YVIEFIAR"};
		peptides=new String[] {"FGGGSVELLK"};

		TDoubleHashSet boundaries=new TDoubleHashSet();
		for (Range range : stripefile.getRanges().keySet()) {
			boundaries.add(range.getStart());
			boundaries.add(range.getStop());
		}
		double[] binArray=boundaries.toArray();
		Arrays.sort(binArray);

		InputStream is=stripefile.getClass().getResourceAsStream("/mouse_20150911_uniprot_sp.fasta");
		ArrayList<FastaEntry> entries=FastaReader.readFasta(is, "mouse_20150911_uniprot_sp.fasta");
		Pair<TDoubleIntHashMap[], ArrayList<String>[]> background=BackgroundGenerator.generateBackground(binArray, entries, PARAMETERS);
		TDoubleIntHashMap[] binCounters=background.x;
		ArrayList<String>[] backgroundProteomes=background.y;
		PrecursorScanMap precursors=new PrecursorScanMap(stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));
		
		PSMScorer scorer=new DotProduct(PARAMETERS.getFragmentTolerance());
		PecanRawScorer pecanScorer=new PecanRawScorer(PARAMETERS.getFragmentTolerance(), new PecanAuxillaryScorer(PARAMETERS, precursors));
		
		// get precursors
		//PecanScorer scorer=new PecanScorer(PARAMETERS.getFragmentTolerance(), PARAMETERS.getPrecursorTolerance(), precursors);
		
		// get stripes
		for (Entry<Range, Float> entry : stripefile.getRanges().entrySet()) {
			Range range=entry.getKey();
			float dutyCycle=entry.getValue();
			int scanAveragingMargin=(int)(PARAMETERS.getMinEluteTime()/dutyCycle/2); // floor
			
			System.out.println("Processing "+range+" ("+scanAveragingMargin+")");
			
			ArrayList<Stripe> stripes=stripefile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, true);
			Collections.sort(stripes);
			
			int index=Arrays.binarySearch(binArray, range.getMiddle());
			index=(-(index+1))-1;
			TDoubleIntHashMap map=binCounters[index];
			double[] keys=map.keys();
			Arrays.sort(keys);
			ArrayList<String> backgroundProteomeArray=backgroundProteomes[index];
			HashSet<String> backgroundProteomeSet=new HashSet<String>(backgroundProteomeArray);
			
			if (peptides!=null) {
				backgroundProteomeArray=new ArrayList<String>(Arrays.asList(peptides));
			}
			
			// first check to see if we need to process this stripe
			boolean hasPeptides=false;
			outer:for (String peptide : backgroundProteomeArray) {
				for (byte charge : charges) {
					double mz=PARAMETERS.getAAConstants().getChargedMass(peptide, charge);
					if (range.contains((float)mz)) {
						hasPeptides=true;
						break outer;
					}
				}
			}
			if (!hasPeptides) {
				continue;
			}

			// prepare executor for background
			ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("SWATH_"+range.getStart()+"to"+range.getStop()+"-%d").setDaemon(true).build();
			LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
			ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 
			//ExecutorService executor=Executors.newFixedThreadPool(cores, threadFactory);

			int backgroundPeptideCount=0;
			int seed=RandomGenerator.randomInt(1);
			ArrayList<Future<HashMap<LibraryEntry, PeptideScoringResult>>> results=new ArrayList<Future<HashMap<LibraryEntry, PeptideScoringResult>>>();
			while (backgroundPeptideCount<2000) {
				for (byte charge : charges) {
					seed=RandomGenerator.randomInt(seed);
					String peptide=backgroundProteomeArray.get((int)(RandomGenerator.floatFromRandomInt(seed)*backgroundProteomeArray.size()));
					double mz=PARAMETERS.getAAConstants().getChargedMass(peptide, charge);

					if (range.contains((float)mz)) {
						String random=PeptideUtils.getDecoy(peptide, backgroundProteomeSet, PARAMETERS);
						FragmentationModel randmodel=new FragmentationModel(random, PARAMETERS.getAAConstants());
						PecanLibraryEntry randentry=randmodel.getPecanSpectrum(charge, keys, map, PARAMETERS);

						ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
						tasks.add(randentry);

						Future<HashMap<LibraryEntry, PeptideScoringResult>> value=executor.submit(new PeptideScoringTask(scorer, tasks, stripes));
						results.add(value);
						
						backgroundPeptideCount++;
					}
				}
			}
			executor.shutdown();
			while (!executor.isTerminated()) {
				System.out.println(workQueue.size()+" background peptides remaining for "+range+"...");
				Thread.sleep(100);
			}
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

			// prepare executor for peptides
			executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 
			
			TDoubleObjectHashMap<TDoubleArrayList> backgroundScoreMap=new TDoubleObjectHashMap<TDoubleArrayList>();
			for (Future<HashMap<LibraryEntry, PeptideScoringResult>> future : results) {
				HashMap<LibraryEntry, PeptideScoringResult> result=future.get();
				for (Entry<LibraryEntry, PeptideScoringResult> resultEntry : result.entrySet()) {
					Pair<double[], double[]> arrays=resultEntry.getValue().getTrace().toArrays();
					double[] x=arrays.x;
					double[] y=arrays.y;
					for (int i=0; i<x.length; i++) {
						TDoubleArrayList list=backgroundScoreMap.get(x[i]);
						if (list==null) {
							list=new TDoubleArrayList();
							backgroundScoreMap.put(x[i], list);
						}
						list.add(y[i]);
					}
				}
			}
			
			final TDoubleObjectHashMap<XYPoint> backgroundScores=new TDoubleObjectHashMap<XYPoint>();
			final ArrayList<XYPoint> meanPlusStdev=new ArrayList<XYPoint>();
			final ArrayList<XYPoint> meanStdev=new ArrayList<XYPoint>();
			final ArrayList<XYPoint> meanMinusStdev=new ArrayList<XYPoint>();
			backgroundScoreMap.forEachEntry(new TDoubleObjectProcedure<TDoubleArrayList>() {
				public boolean execute(double arg0, TDoubleArrayList arg1) {
					double[] values=arg1.toArray();
					double m=General.mean(values);
					double s=General.stdev(values);
					backgroundScores.put(arg0, new XYPoint(m, s));

					meanPlusStdev.add(new XYPoint(arg0, m+s));
					meanStdev.add(new XYPoint(arg0, m));
					meanMinusStdev.add(new XYPoint(arg0, m-s));
					return true;
				};
			});
			//Charter.launchChart("RT ("+range+" M/Z)", "Background Score", new XYTrace(meanPlusStdev, GraphType.line, "M+S"), new XYTrace(meanStdev, GraphType.line, "M"), new XYTrace(meanMinusStdev, GraphType.line, "M-S"));

			results.clear();
			for (String peptide : backgroundProteomeArray) {
				for (byte charge : charges) {
					double mz=PARAMETERS.getAAConstants().getChargedMass(peptide, charge);
					if (range.contains((float)mz)) {
						FragmentationModel model=new FragmentationModel(peptide, PARAMETERS.getAAConstants());
						PecanLibraryEntry pecanEntry=model.getPecanSpectrum(charge, keys, map, PARAMETERS);

						FragmentationModel revmodel=new FragmentationModel(PeptideUtils.getSmartDecoy(peptide, backgroundProteomeSet, PARAMETERS), PARAMETERS.getAAConstants());
						PecanLibraryEntry reventry=revmodel.getPecanSpectrum(charge, keys, map, PARAMETERS);

						ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
						tasks.add(pecanEntry);
						tasks.add(reventry);

						Future<HashMap<LibraryEntry, PeptideScoringResult>> value=executor.submit(new PecanScoringTask(pecanScorer, tasks, stripes, backgroundScores, scanAveragingMargin));
						results.add(value);
					}
				}
			}
			executor.shutdown();
			while (!executor.isTerminated()) {
				System.out.println(workQueue.size()+" peptides remaining for "+range+"...");
				Thread.sleep(100);
			}
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

			ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
			for (Future<HashMap<LibraryEntry, PeptideScoringResult>> future : results) {
				HashMap<LibraryEntry, PeptideScoringResult> result=future.get();
				for (Entry<LibraryEntry, PeptideScoringResult> resultEntry : result.entrySet()) {
					LibraryEntry peptide=resultEntry.getKey();
					PeptideScoringResult peptideResult=resultEntry.getValue();
					
					int rank=1;
					for (Pair<ScoredObject<Stripe>, float[]> goodStripe : peptideResult.getGoodStripes()) {
						float primaryScore=goodStripe.x.x;
						Stripe stripe=goodStripe.x.y;
						float[] auxScores=goodStripe.y;
						
						System.out.print(peptide.getPeptideModSeq()+"\t"+rank+"\t"+primaryScore+"\t"+stripe.getScanStartTime());
						for (float s : auxScores) {
							System.out.print("\t"+s);
						}
						System.out.println();
						rank++;
						if (rank>3) break;
					}
					traces.add(peptideResult.getTrace());
				}
			}
			//Charter.launchChart("RT ("+range+" M/Z)", "Score", traces.toArray(new XYTrace[traces.size()]));
		}
	}

}
