package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ExecutionException;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PecanLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import gnu.trove.map.hash.TDoubleIntHashMap;
import gnu.trove.set.hash.TDoubleHashSet;
import junit.framework.TestCase;

public class PeptideScoringTaskTest extends TestCase {
	private static final SearchParameters PARAMETERS=new SearchParameters(FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"));
	
	public void testPeptideScoringTask() throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {

		File f=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/82593_lv_mcx_DIA_5mz_400to525.dia");
		StripeFile stripefile=new StripeFile();
		stripefile.openFile(f);

		String peptide="FGGGSVELLK";
		byte charge=(byte)2;

		for (Range range : stripefile.getRanges().keySet()) {
			// first check to see if we need to process this stripe
			boolean hasPeptides=false;
			double mz=MassConstants.getChargedMass(peptide, charge);
			if (!range.contains((float)mz)) {
				continue;
			}

			TDoubleHashSet boundaries=new TDoubleHashSet();
			boundaries.add(range.getStart());
			boundaries.add(range.getStop());
			double[] binArray=boundaries.toArray();
			Arrays.sort(binArray);

			InputStream is=stripefile.getClass().getResourceAsStream("/mouse_20150911_uniprot_sp.fasta");
			ArrayList<FastaEntry> entries=FastaReader.readFasta(is, "mouse_20150911_uniprot_sp.fasta");
			Pair<TDoubleIntHashMap[], ArrayList<String>[]> background=BackgroundGenerator.generateBackground(binArray, entries, PARAMETERS);
			TDoubleIntHashMap[] binCounters=background.x;
			ArrayList<String>[] backgroundProteomes=background.y;
			PSMScorer scorer=new DotProduct(PARAMETERS.getFragmentTolerance());
			
			ArrayList<Stripe> stripes=stripefile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, true);
			int index=Arrays.binarySearch(binArray, range.getMiddle());
			index=(-(index+1))-1;
			TDoubleIntHashMap map=binCounters[index];
			double[] keys=map.keys();
			Arrays.sort(keys);
			ArrayList<String> backgroundProteomeArray=backgroundProteomes[index];
			HashSet<String> backgroundProteomeSet=new HashSet<String>(backgroundProteomeArray);
			
			FragmentationModel model=new FragmentationModel(peptide);
			PecanLibraryEntry entry=model.getPecanSpectrum(charge, keys, map, PARAMETERS);

			FragmentationModel revmodel=new FragmentationModel(PeptideUtils.getSmartDecoy(peptide, backgroundProteomeSet, PARAMETERS));
			PecanLibraryEntry reventry=revmodel.getPecanSpectrum(charge, keys, map, PARAMETERS);

			ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
			tasks.add(entry);

			PeptideScoringTask task=new PeptideScoringTask(scorer, tasks, stripes);
			
			HashMap<LibraryEntry, PeptideScoringResult> result=task.process();
		}
	}

}
