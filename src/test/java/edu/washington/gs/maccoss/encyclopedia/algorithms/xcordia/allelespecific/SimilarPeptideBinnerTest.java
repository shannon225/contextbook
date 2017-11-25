package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import junit.framework.TestCase;

public class SimilarPeptideBinnerTest extends TestCase {
	public static void main(String[] args) throws Exception {
		File peffFile=new File("/Users/searleb/Documents/school/xcordia_manuscript/amyloid_protein.peff");
		InputStream is=new FileInputStream(peffFile);
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(new BufferedReader(new InputStreamReader(is)), peffFile.getName(), null, true);
		PecanSearchParameters parameters=PecanParameterParser.getDefaultParametersObject();
		
		HashSet<FastaPeptideEntry> targets=new HashSet<>();
		for (FastaEntryInterface protein : entries) {
			ArrayList<String> peptides=parameters.getEnzyme().digestProtein(protein, 6, 100, 0, parameters.getAAConstants());
			for (String peptide : peptides) {
				FastaPeptideEntry pe=protein.getSubEntry(peptide);
				targets.add(pe);
			}
		}
		
		SimilarPeptideBinner binner=new SimilarPeptideBinner();
		ArrayList<ArrayList<FastaPeptideEntry>> bins=binner.binPeptides(targets);
		System.out.println(targets.size()+"/"+bins.size());
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
