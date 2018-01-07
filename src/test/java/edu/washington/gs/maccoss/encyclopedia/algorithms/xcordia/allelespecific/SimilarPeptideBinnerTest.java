package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.XCorDIA;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.BackgroundFrequencyInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.UnitBackgroundFrequencyCalculator;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideDatabase;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import junit.framework.TestCase;

public class SimilarPeptideBinnerTest extends TestCase {
	public static void main(String[] args) throws Exception {

		HashMap<String, String> defaults=PecanParameterParser.getDefaultParameters();
		defaults.put("-localizationModification", "Phosphorylation");
		defaults.put("-scoringBreadthType", "uncal20");
		PecanSearchParameters parameters=PecanParameterParser.parseParameters(defaults);
		
		System.out.println("Reading raw file...");
		File diaFile=new File("/Users/searleb/Documents/school/xcordia_manuscript/demux/20141121_3_4_DIA_1.dia");
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
		
		UnitBackgroundFrequencyCalculator unitBackgroundFrequencyCalculator=new UnitBackgroundFrequencyCalculator(0.01f);
		BackgroundFrequencyInterface background=unitBackgroundFrequencyCalculator;
		//background=BackgroundFrequencyCalculator.generateBackground(stripefile);
		
		PhosphoLocalizer localizer=new PhosphoLocalizer(stripefile, PeptideModification.polymorphism, background, parameters);
		
		ArrayList<Range> ranges=new ArrayList<>(stripefile.getRanges().keySet());
		Collections.sort(ranges);

		System.out.println("Reading peff fasta file...");
		File peffFile=new File("/Users/searleb/Documents/school/xcordia_manuscript/amyloid_protein.peff");
		InputStream is=new FileInputStream(peffFile);
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(new BufferedReader(new InputStreamReader(is)), peffFile.getName(), null, true);
		
		PeptideDatabase targets=new PeptideDatabase();
		for (FastaEntryInterface protein : entries) {
			ArrayList<String> peptides=parameters.getEnzyme().digestProtein(protein, 6, 100, 0, parameters.getAAConstants());
			for (String peptide : peptides) {
				FastaPeptideEntry pe=protein.getSubEntry(peptide);
				targets.add(pe);
			}
		}
		
		System.out.println("Total unique peptides: "+targets.size());
		for (Range range : ranges) {
			HashSet<FastaPeptideEntry> peptides=XCorDIA.getPeptidesInRange(parameters, targets, range);
			SimilarPeptideBinner binner=new SimilarPeptideBinner();
			ArrayList<ArrayList<FastaPeptideEntry>> bins=binner.binPeptides(peptides);
			int[] counts=new int[6];
			for (ArrayList<FastaPeptideEntry> bin : bins) {
				int index=Math.min(counts.length-1, bin.size());
				if (bin.size()<counts.length) {
					counts[index]++;
				}
				//PhosphoLocalizationData actuallyPhosphoData=localizer.extractPhosphoFormsFromStripes(peptideModSeq, precursorMz, precursorCharge, permutations, retentionTime, stripes, true);
			}
			System.out.println(range.toString()+"\t"+peptides.size()+"\t"+bins.size()+"\t"+(peptides.size()-bins.size())+"\t"+General.toString(counts));
		}
	}
	
	public void testBinner() {
		String annotation=">nxp:NX_P0DJI8-1 \\DbUniqueId=NX_P0DJI8-1 \\PName=Serum amyloid A-1 protein isoform Iso 1 \\GName=SAA1 \\NcbiTaxId=9606 \\TaxName=Homo Sapiens \\Length=122 \\SV=95 \\EV=159 \\PE=1 \\ModResPsi=(101|MOD:00316|N4,N4-dimethyl-L-asparagine) \\VariantSimple=(15|S)(70|A)(75|V)(78|N)(86|L)(90|D) \\Processed=(1|18|signal peptide)(19|94|mature protein)(19|122|mature protein)(20|120|mature protein)(20|121|mature protein)(20|122|mature protein)(21|122|mature protein)(22|119|mature protein)(95|122|maturation peptide)";
		String simpleAnnotation=">nxp:NX_P0DJI8-1";
		String sequence="MKLLTGLVFCSLVLGVSSRSFFSFLGEAFDGARDMWRAYSDMREANYIGSDKYFHARGNYDAAKRGPGGVWAAEAISDARENIQRFFGHGAEDSLADQAANEWGRSGKDPNHFRPAGLPEKY";
		FastaEntry simpleEntry=new FastaEntry("source", simpleAnnotation, sequence);
		ExtendedFastaEntry entry=new ExtendedFastaEntry("source", annotation, sequence);
		
		PecanSearchParameters parameters=PecanParameterParser.getDefaultParametersObject();
		ArrayList<String> simplePeptides=parameters.getEnzyme().digestProtein(simpleEntry, 6, 100, 0, parameters.getAAConstants());
		ArrayList<String> peptides=parameters.getEnzyme().digestProtein(entry, 6, 100, 0, parameters.getAAConstants());

		HashSet<FastaPeptideEntry> targets=new HashSet<>();
		for (String peptide : peptides) {
			FastaPeptideEntry pe=entry.getSubEntry(peptide);
			targets.add(pe);
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
