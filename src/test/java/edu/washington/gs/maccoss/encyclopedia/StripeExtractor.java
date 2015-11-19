package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
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
import edu.washington.gs.maccoss.encyclopedia.algorithms.PecanRawScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringTask;
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
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TDoubleIntHashMap;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.procedure.TDoubleObjectProcedure;
import gnu.trove.set.hash.TDoubleHashSet;

public class StripeExtractor {
	private static final SearchParameters PARAMETERS=new SearchParameters(FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"));

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
		String[] peptides=new String[] {"FGGGSVELLK"};
		//String[] peptides=new String[] {"ADVDAATLAR", "APIQWEER", "ATNLTVSAVR", "AVDSLVPIGR", "AYIDSTDSR", "DGLTDVYNK", "DGPGFYTTR", "DTPGFIVNR", "DVLSNLIPK", "DYPLIPVGK", "EAFLLFDR", "EALISQLTR",
		//		"FGGGSVELLK", "GGASDALLYR", "GIWHNYDK", "GLNEAAIVNK", "GPLVQGVDSR", "GSGIQWDLR", "GTAVVNGEFK", "IGGIGTVPVGR", "IGYPAPNFK", "ILQEGVDPK", "IRCDIANVK", "ITQSNAILR", "IVVHAGGVIR",
		//		"KFVADGIFK", "LGFMSAFVK", "LLEAASVSSK", "LLFEELVR", "LNVLANVIR", "LQGDLVTIR", "LTLSALIDGK", "LTLSALVDGK", "LVNMLDAVR", "MFASFPTTK", "NDAGYSEPR", "NTYYASIAK", "QGVLTLEIR", "QIFLGGVDR",
		//		"RAEVLDSTK", "RGDFIPGLR", "RLTDADAMK", "RLVVQQAGK", "RWEVAALR", "SFLPLLRR", "SLHTLFGDK", "SVQAAMEKR", "TAYVGENVR", "TDLTAVPASR", "TIPWLENR", "TLEDILFR", "TQLVSNLKK", "TQVQSVIDK",
		//		"TSGGAGGLGSLR", "VAAENQYGR", "VDFDDIHR", "VGPANPSLQK", "VLSIGDGIAR", "VTLVSAAPEK", "VVDDELATR", "VVFIFGPDK", "WPLYLSTK", "YDHLGDSPK", "YVDMSAKSK", "YVIEFIAR"};

		TDoubleHashSet boundaries=new TDoubleHashSet();
		for (Range range : stripefile.getRanges()) {
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
		PSMScorer scorer=new DotProduct(PARAMETERS.getFragmentTolerance());
		
		// get precursors
		PrecursorScanMap precursors=new PrecursorScanMap(stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));
		//PecanScorer scorer=new PecanScorer(PARAMETERS.getFragmentTolerance(), PARAMETERS.getPrecursorTolerance(), precursors);
		
		// get stripes
		for (Range range : stripefile.getRanges()) {
			System.out.println("Processing "+range);
			// first check to see if we need to process this stripe
			boolean hasPeptides=false;
			outer:for (String peptide : peptides) {
				for (byte charge : charges) {
					double mz=MassConstants.getChargedMass(peptide, charge);
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
			
			ArrayList<Stripe> stripes=stripefile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, true);
			int index=Arrays.binarySearch(binArray, range.getMiddle());
			index=(-(index+1))-1;
			TDoubleIntHashMap map=binCounters[index];
			double[] keys=map.keys();
			Arrays.sort(keys);
			ArrayList<String> backgroundProteomeArray=backgroundProteomes[index];
			HashSet<String> backgroundProteomeSet=new HashSet<String>(backgroundProteomeArray);

			int backgroundPeptideCount=0;
			int seed=RandomGenerator.randomInt(1);
			ArrayList<Future<HashMap<LibraryEntry, XYTrace>>> results=new ArrayList<Future<HashMap<LibraryEntry, XYTrace>>>();
			while (backgroundPeptideCount<20) {
				for (byte charge : charges) {
					seed=RandomGenerator.randomInt(seed);
					String peptide=backgroundProteomeArray.get((int)(RandomGenerator.floatFromRandomInt(seed)*backgroundProteomeArray.size()));
					double mz=MassConstants.getChargedMass(peptide, charge);

					if (range.contains((float)mz)) {
						String random=PeptideUtils.getDecoy(peptide, PARAMETERS.getEnzyme(), backgroundProteomeSet);
						FragmentationModel randmodel=new FragmentationModel(random);
						PecanLibraryEntry randentry=randmodel.getPecanSpectrum(charge, keys, map, PARAMETERS);

						ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
						tasks.add(randentry);

						Future<HashMap<LibraryEntry, XYTrace>> value=executor.submit(new PeptideScoringTask(scorer, tasks, stripes));
						results.add(value);
						
						backgroundPeptideCount++;
					}
				}
			}
			executor.shutdown();
			while (!executor.isTerminated()) {
				System.out.println(workQueue.size()+" peptides remaining...");
				Thread.sleep(50);
			}
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

			// prepare executor for peptides
			executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 
			
			TDoubleObjectHashMap<TDoubleArrayList> backgroundScoreMap=new TDoubleObjectHashMap<TDoubleArrayList>();
			for (Future<HashMap<LibraryEntry, XYTrace>> future : results) {
				HashMap<LibraryEntry, XYTrace> result=future.get();
				for (Entry<LibraryEntry, XYTrace> entry : result.entrySet()) {
					Pair<double[], double[]> arrays=entry.getValue().toArrays();
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
			Charter.launchChart("RT ("+range+" M/Z)", "Background Score", new XYTrace(meanPlusStdev, GraphType.line, "M+S"), new XYTrace(meanStdev, GraphType.line, "M"), new XYTrace(meanMinusStdev, GraphType.line, "M-S"));

			scorer=new PecanRawScorer(PARAMETERS.getFragmentTolerance());
			results.clear();
			for (String peptide : peptides) {
				for (byte charge : charges) {
					double mz=MassConstants.getChargedMass(peptide, charge);
					if (range.contains((float)mz)) {
						FragmentationModel model=new FragmentationModel(peptide);
						PecanLibraryEntry entry=model.getPecanSpectrum(charge, keys, map, PARAMETERS);

						FragmentationModel revmodel=new FragmentationModel(PeptideUtils.getSmartDecoy(peptide, backgroundProteomeSet, PARAMETERS));
						PecanLibraryEntry reventry=revmodel.getPecanSpectrum(charge, keys, map, PARAMETERS);

						ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
						tasks.add(entry);
						tasks.add(reventry);

						Future<HashMap<LibraryEntry, XYTrace>> value=executor.submit(new PeptideScoringTask(scorer, tasks, stripes, backgroundScores));
						results.add(value);
					}
				}
			}
			executor.shutdown();
			while (!executor.isTerminated()) {
				executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
			}

			ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
			for (Future<HashMap<LibraryEntry, XYTrace>> future : results) {
				HashMap<LibraryEntry, XYTrace> result=future.get();
				for (Entry<LibraryEntry, XYTrace> entry : result.entrySet()) {
					traces.add(entry.getValue());
				}
			}
			Charter.launchChart("RT ("+range+" M/Z)", "Score", traces.toArray(new XYTrace[traces.size()]));
		}
	}

}
