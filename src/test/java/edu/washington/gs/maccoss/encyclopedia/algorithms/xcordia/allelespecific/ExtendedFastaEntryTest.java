package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import junit.framework.TestCase;
import uk.ac.ebi.pride.utilities.pridemod.ModReader;
import uk.ac.ebi.pride.utilities.pridemod.model.PTM;
import uk.ac.ebi.pride.utilities.pridemod.model.Specificity;

public class ExtendedFastaEntryTest extends TestCase {
	public static void main(String[] args) {
		ModReader modReader = ModReader.getInstance();
		DecimalFormat myFormatter = new DecimalFormat("00000");
		for (int i = 1; i < 100; i++) {
			String accession="MOD:" + myFormatter.format(i);
			PTM ptm = modReader.getPTMbyAccession(accession);
			if (ptm==null) {
				System.out.println("No mod for "+accession);
			} else {
				System.out.print(accession+", "+ptm.getName() + ", " + ptm.getMonoDeltaMass()+": ");
				for (Specificity spec : ptm.getSpecificityCollection()) {
					System.out.print(spec.getName());
				}
				System.out.println();
			}
		}
	}

	/* peff format
	 * 
	 * VariantSimple=(position|newAminoAcid|optionalTag)
	 * (ambiguity codes such as J or X are permitted)
	 * * = stop codon
	 * only for single amino acid substitutions
	 * 
	 * e.g. VariantSimple=(223|A)(225|C|dbSNP)
	 * (startPosition|endPosition|newSequence|optionalTag)
	 * only the cases that can not be encoded using VariantSimple
	 * this is for insertions or deletions
	 * 
	 * e.g.: VariantComplex=(4|4|)(32|32|LM)(100|101|P)(100|102|KPA)(100|100||10kexomes)
	 */
	
	public void testParsing() {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		String annotation=">nxp:NX_P02766-1 \\DbUniqueId=NX_P02766-1 \\PName=Transthyretin isoform Iso 1 \\GName=TTR \\NcbiTaxId=9606 \\TaxName=Homo Sapiens \\Length=147 \\SV=266 \\EV=656 \\PE=1 \\ModResPsi=(62|MOD:00041|L-gamma-carboxyglutamic acid)(69|MOD:00047|O-phospho-L-threonine)(72|MOD:00046|O-phospho-L-serine) \\ModRes=(118||N-linked (GlcNAc...)) "
				+"\\VariantSimple=(2|G)(5|H)(5|C)(8|H)(9|V)(9|F)(13|E)(13|R)(18|F)(18|A)(19|D)(23|M)(26|D)(26|S)(30|R)(32|P)(33|T)(33|I)(38|E)(38|G)(40|I)(41|Q)(42|D)(43|N)(44|S)(45|T)(46|V)(47|S)(48|M)(50|L)(50|G)(50|A)(50|M)(51|R)(53|C)(53|V)(53|L)(53|I)(54|T)(55|N)(55|M)(56|P)(58|A)(58|V)(61|L)(62|D)(62|G)(64|L)(64|S)(65|S)(65|D)(65|T)(67|E)(67|V)(67|R)(67|A)(69|I)(69|A)(70|N)(70|G)(70|R)(70|I)(72|P)(72|Y)(73|E)(74|G)(74|K)(74|S)(75|Q)(75|P)(76|Y)(78|H)(78|R)(79|I)(79|A)(79|K)(80|A)(81|G)(81|K)(82|D)(84|L)(88|L)(89|H)(90|N)(91|A)(93|V)(94|H)(95|S)(97|Y)(98|F)(99|*)(100|E)(104|N)(104|T)(104|L)(104|S)(105|P)(109|K)(109|D)(109|Q)(110|N)(111|S)(112|*)(114|A)(117|G)(117|S)(119|A)(121|D)(121|S)(121|A)(122|R)(123|H)(124|C)(124|H)(125|H)(125|*)(126|N)(126|I)(127|M)(127|V)(129|T)(131|M)(134|H)(134|C)(136|V)(136|H)(136|S)(139|M)(140|V)(140|S)(142|A)(142|I)(144|S)(145|S)(52|E)(80|I)(84|S)(89|C)(90|Q)(101|S)(101|T)(103|A)(118|K)(119|N)(120|P)(123|S)(123|C)(127|T)(129|V)(133|S)(135|F)(136|C)(137|T)(146|R)(147|D)(147|*) "
				+"\\VariantComplex=(3|3|SH)(4|4|)(32|32|LM)(142|142|) \\Processed=(1|20|signal peptide)(21|147|mature protein)\n";
		String sequence="MASHRLLLLCLAGLVFVSEAGPTGTGESKCPLMVKVLDAVRGSPAINVAVHVFRKAADDTWEPFASGKTSESGELHGLTTEEEFVEGIYKVEIDTKSYWKALGISPFHEHAEVVFTANDSGPRRYTIAALLSPYSYSTTAVVTNPKE";
		ExtendedFastaEntry extendedEntry=new ExtendedFastaEntry("testExtendedFastaEntry", annotation, sequence, parameters);
		ArrayList<AlleleVariant> variants=extendedEntry.getPotentialVariants();
		
		assertEquals(158,variants.size());
		assertEquals(extendedEntry.getSequence(),sequence);
		
		for (int index=0; index<variants.size(); index++) {
			AlleleVariant variant=variants.get(index);
			assertEquals(sequence.substring(variant.getStartSite()-1, variant.getStopSite()),variant.getOriginalSequence());
			String newSeq = sequence.substring(0,variant.getStartSite()-1)+variant.getNewSequence()+sequence.substring(variant.getStopSite());
			assertEquals(newSeq.length(),sequence.length()-variant.getOriginalSequence().length() +variant.getNewSequence().length());
		}
	}
	public void testParsingMods() {
		String annotation=">nxp:NX_Q96B36-1 \\DbUniqueId=NX_Q96B36-1 \\PName=Proline-rich AKT1 substrate 1 isoform Iso 1 \\GName=AKT1S1 \\NcbiTaxId=9606 \\TaxName=Homo Sapiens \\Length=256 \\SV=135 \\EV=428 \\PE=1 "
				+ "\\ModResPsi=(51|MOD:00078|omega-N-methyl-L-arginine)(73|MOD:00047|O-phospho-L-threonine)(90|MOD:00047|O-phospho-L-threonine)(97|MOD:00047|O-phospho-L-threonine)(116|MOD:00046|O-phospho-L-serine)(187|MOD:00046|O-phospho-L-serine)(198|MOD:00047|O-phospho-L-threonine)(247|MOD:00046|O-phospho-L-serine)(88|MOD:00046|O-phospho-L-serine)(92|MOD:00046|O-phospho-L-serine)(183|MOD:00046|O-phospho-L-serine)(202|MOD:00046|O-phospho-L-serine)(203|MOD:00046|O-phospho-L-serine)(211|MOD:00046|O-phospho-L-serine)(212|MOD:00046|O-phospho-L-serine)(246|MOD:00047|O-phospho-L-threonine) "
				+ "\\VariantSimple=(1|V)(3|L)(5|H)(12|T)(13|M)(19|C)(21|W)(21|Q)(23|Q)(26|M)(27|D)(33|T)(40|H)(41|L)(42|S)(43|S)(46|H)(46|C)(47|P)(51|Q)(51|L)(55|V)(59|H)(60|C)(60|H)(61|Y)(63|R)(68|G)(76|Q)(76|W)(78|A)(79|V)(80|S)(82|Q)(86|L)(87|S)(88|C)(88|R)(93|T)(95|W)(95|Q)(95|P)(97|A)(99|V)(103|D)(106|H)(107|K)(114|G)(115|I)(122|T)(122|V)(129|I)(130|A)(135|T)(136|I)(141|H)(141|S)(142|R)(142|H)(142|S)(145|K)(147|E)(149|K)(151|I)(161|S)(162|T)(163|D)(166|I)(168|L)(170|H)(172|V)(174|T)(182|E)(185|S)(186|L)(188|L)(190|F)(190|I)(192|S)(197|T)(198|I)(201|W)(201|Q)(204|N)(207|D)(209|R)(209|L)(214|Y)(215|P)(223|H)(233|A)(233|S)(234|H)(239|V)(242|L)(245|K)(117|R)(212|L)(249|L)(252|Q)(7|K)(51|*)(53|S)(58|V)(64|N)(104|K)(158|*)(163|A)(199|*)(200|V) \\VariantComplex=(39|39|PP)(206|206|)(114|115|) \\Processed=(1|256|mature protein)";
		String sequence="MASGRPEELWEAVVGAAERFRARTGTELVLLTAAPPPPPRPGPCAYAAHGRGALAEAARR" + 
				"CLHDIALAHRAATAARPPAPPPAPQPPSPTPSPPRPTLAREDNEEDEDEPTETETSGEQL" + 
				"GISDNGGLFVMDEDATLQDLPPFCESDPESTDDGSLSEETPAGPPTCSVPPASALPTQQY" + 
				"AKSLPVSVPVWGFKEKRTEARSSDEENGPPSSPDLDRIAASMRALVLREAEDTQVFGDLP" + 
				"RPRLNTSDFQKLKRKY";

		// no mods (so use variants)
		HashMap<String, String> map=PecanParameterParser.getDefaultParameters();
		SearchParameters parameters=PecanParameterParser.parseParameters(map);
		ExtendedFastaEntry extendedEntry=new ExtendedFastaEntry("testExtendedFastaEntry", annotation, sequence, parameters);
		ArrayList<AlleleVariant> variants=extendedEntry.getPotentialVariants();
		
		assertEquals(113,variants.size());
	
		// 1234567890 1234567890 1234567890 1234567890 1234567890 1234567890
		// MASGRPEELW EAVVGAAERF RARTGTELVL LTAAPPPPPR PGPCAYAAHG RGALAEAARR 60
		// CLHDIALAHR AATAARPPAP PPAPQPPSPT PSPPRPTLAR EDNEEDEDEP TETETSGEQL 120
		// GISDNGGLFV MDEDATLQDL PPFCESDPES TDDGSLSEET PAGPPTCSVP PASALPTQQY 180
		// AKSLPVSVPV WGFKEKRTEA RSSDEENGPP SSPDLDRIAA SMRALVLREA EDTQVFGDLP 240
		// RPRLNTSDFQ KLKRKY 

		// phosphorylation (no variants)
		map.put("-variable", "S=79.966331,T=79.966331,Y=79.966331");
		parameters=PecanParameterParser.parseParameters(map);
		extendedEntry=new ExtendedFastaEntry("testExtendedFastaEntry", annotation, sequence, parameters);
		variants=extendedEntry.getPotentialVariants();
		assertEquals(16-1,variants.size());
		
		// methylation (no variants)
		map.put("-variable", "R=14.015650,K=14.015650,H=14.015650");
		parameters=PecanParameterParser.parseParameters(map);
		extendedEntry=new ExtendedFastaEntry("testExtendedFastaEntry", annotation, sequence, parameters);
		variants=extendedEntry.getPotentialVariants();
		assertEquals(16-15,variants.size());
	}
}
