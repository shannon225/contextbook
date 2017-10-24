package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific;

import java.util.ArrayList;
import junit.framework.TestCase;

public class ExtendedFastaEntryTest extends TestCase {
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
		String annotation=">nxp:NX_P02766-1 \\DbUniqueId=NX_P02766-1 \\PName=Transthyretin isoform Iso 1 \\GName=TTR \\NcbiTaxId=9606 \\TaxName=Homo Sapiens \\Length=147 \\SV=266 \\EV=656 \\PE=1 \\ModResPsi=(62|MOD:00041|L-gamma-carboxyglutamic acid)(69|MOD:00047|O-phospho-L-threonine)(72|MOD:00046|O-phospho-L-serine) \\ModRes=(118||N-linked (GlcNAc...)) "
				+"\\VariantSimple=(2|G)(5|H)(5|C)(8|H)(9|V)(9|F)(13|E)(13|R)(18|F)(18|A)(19|D)(23|M)(26|D)(26|S)(30|R)(32|P)(33|T)(33|I)(38|E)(38|G)(40|I)(41|Q)(42|D)(43|N)(44|S)(45|T)(46|V)(47|S)(48|M)(50|L)(50|G)(50|A)(50|M)(51|R)(53|C)(53|V)(53|L)(53|I)(54|T)(55|N)(55|M)(56|P)(58|A)(58|V)(61|L)(62|D)(62|G)(64|L)(64|S)(65|S)(65|D)(65|T)(67|E)(67|V)(67|R)(67|A)(69|I)(69|A)(70|N)(70|G)(70|R)(70|I)(72|P)(72|Y)(73|E)(74|G)(74|K)(74|S)(75|Q)(75|P)(76|Y)(78|H)(78|R)(79|I)(79|A)(79|K)(80|A)(81|G)(81|K)(82|D)(84|L)(88|L)(89|H)(90|N)(91|A)(93|V)(94|H)(95|S)(97|Y)(98|F)(99|*)(100|E)(104|N)(104|T)(104|L)(104|S)(105|P)(109|K)(109|D)(109|Q)(110|N)(111|S)(112|*)(114|A)(117|G)(117|S)(119|A)(121|D)(121|S)(121|A)(122|R)(123|H)(124|C)(124|H)(125|H)(125|*)(126|N)(126|I)(127|M)(127|V)(129|T)(131|M)(134|H)(134|C)(136|V)(136|H)(136|S)(139|M)(140|V)(140|S)(142|A)(142|I)(144|S)(145|S)(52|E)(80|I)(84|S)(89|C)(90|Q)(101|S)(101|T)(103|A)(118|K)(119|N)(120|P)(123|S)(123|C)(127|T)(129|V)(133|S)(135|F)(136|C)(137|T)(146|R)(147|D)(147|*) "
				+"\\VariantComplex=(3|3|SH)(4|4|)(32|32|LM)(142|142|) \\Processed=(1|20|signal peptide)(21|147|mature protein)\n";
		String sequence="MASHRLLLLCLAGLVFVSEAGPTGTGESKCPLMVKVLDAVRGSPAINVAVHVFRKAADDTWEPFASGKTSESGELHGLTTEEEFVEGIYKVEIDTKSYWKALGISPFHEHAEVVFTANDSGPRRYTIAALLSPYSYSTTAVVTNPKE";
		System.out.println("PEFF format variant annotation parsing test");
		ExtendedFastaEntry extendedEntry=new ExtendedFastaEntry("testExtendedFastaEntry", annotation, sequence);
		ArrayList<AlleleVariant> variants=extendedEntry.getPotentialVariant();
		
		assertEquals(160,variants.size());
		
		String entrySequence=extendedEntry.getSequence();
		for (int index=0; index<variants.size(); index++) {
			AlleleVariant variant=variants.get(index);
			assertEquals(entrySequence.substring(variant.getStartSite()-1, variant.getStopSite()),variant.getOriginalSequence());
			
		}
	}
	
}
