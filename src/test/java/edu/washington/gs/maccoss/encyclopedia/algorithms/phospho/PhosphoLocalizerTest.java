package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.awt.Color;
import java.awt.Dimension;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.ChromatogramExtractor;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.map.hash.TFloatFloatHashMap;
import junit.framework.TestCase;

public class PhosphoLocalizerTest extends TestCase {
	public static final Color[] colors=new Color[] {Color.red, Color.blue, Color.green, Color.cyan, Color.magenta}; 
	public static void main(String[] args) throws Exception {
		File libraryFile=new File("/Users/searleb/Documents/school/projects/may_asms/hela_phospho/VillenJ_Exactive_HumanPhosphoproteome.elib");
		File diaFile=new File("/Users/searleb/Documents/school/projects/may_asms/hela_phospho/110515_bcs_hela_phospho_starved_20mz_500_900.dia");
		//diaFile=new File("/Users/searleb/Documents/school/projects/mzml/q06048_rl_MCF7_IMAC_GpX_3.dia");

		//libraryFile=new File("/Users/searleb/Documents/school/projects/test.elib");
		//diaFile=new File("/Users/searleb/Documents/school/projects/may_asms/phospho/22may2016_mcf7_dia_reserveA_01.mzML");

		//libraryFile=new File("/Users/searleb/Documents/school/projects/test.elib");
		//diaFile=new File("/Users/searleb/Documents/school/projects/may_asms/phospho/22may2016_mcf7_dia_reserveA_01.mzML");


		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		
		LibraryFile library=new LibraryFile();
		library.openFile(libraryFile);
		
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
		
		PhosphoLocalizer localizer=new PhosphoLocalizer(stripefile, library, parameters);

		PSMData psmdata=getPeptide(10);
		//PSMData psmdata=getPeptide(0);
		
		ArrayList<Spectrum> precursors=new ArrayList<Spectrum>();

		float duration=1.5f*60f; // search for 5 minutes
		for (PrecursorScan stripe : stripefile.getPrecursors(psmdata.getRetentionTime()-duration, psmdata.getRetentionTime()+duration)) {
			precursors.add(stripe);
		}
		Charter.launchChart("Retention Time", "Intensity", false, new Dimension(800, 250), ChromatogramExtractor.extractPrecursorChromatograms(parameters.getPrecursorTolerance(), psmdata.getPrecursorMZ(), psmdata.getPrecursorCharge(), precursors));
		//ArrayList<String> permutations=new ArrayList<String>();
		//permutations.add("DKRPLS[+79.96633]GPDVGTPQPAGLASGAK");
		//permutations.add("DKRPLSGPDVGTPQPAGLAS[+79.96633]GAK");
		//ArrayList<Stripe> stripes=stripefile.getStripes(psmdata.getPrecursorMZ(), 0, Float.MAX_VALUE, false);
		//PhosphoLocalizationData multiple=localizer.extractPhosphoForms(psmdata.getPrecursorMZ(), psmdata.getPrecursorCharge(), permutations, psmdata.getRetentionTime(), stripes);
		ArrayList<Stripe> stripes=stripefile.getStripes(psmdata.getPrecursorMZ(), 0, Float.MAX_VALUE, false);
		HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> allVsUniqueList=localizer.runDIAPhosphoLocalization(psmdata, stripes).get().getScoreTraces();
		//HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> allVsUniqueList=multiple.getTraces();
		
		ArrayList<Color> shades=new ArrayList<Color>(Arrays.asList(colors));
		ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
		for (Entry<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> entry : allVsUniqueList.entrySet()) {
			String seq=entry.getKey();
			Pair<TFloatFloatHashMap, TFloatFloatHashMap> pair=entry.getValue();
			Color color=shades.remove(0);
			//traces.add(new XYTrace(pair.x, GraphType.line, "ALL_"+seq, color, 5.0f));
			traces.add(new XYTrace(pair.y, GraphType.line, "UNI_"+seq, color, 3.0f));
		}
		
		Charter.launchChart("Retention Time (Site Specific)", "Score", false, new Dimension(800, 250), traces.toArray(new XYTrace[traces.size()]));

		shades=new ArrayList<Color>(Arrays.asList(colors));
		traces=new ArrayList<XYTrace>();
		for (Entry<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> entry : allVsUniqueList.entrySet()) {
			String seq=entry.getKey();
			Pair<TFloatFloatHashMap, TFloatFloatHashMap> pair=entry.getValue();
			Color color=shades.remove(0);
			traces.add(new XYTrace(pair.x, GraphType.line, "ALL_"+seq, color, 3.0f));
			//traces.add(new XYTrace(pair.y, GraphType.line, "UNI_"+seq, color, 5.0f));
		}
		
		Charter.launchChart("Retention Time (All Ions)", "Score", false, new Dimension(800, 250), traces.toArray(new XYTrace[traces.size()]));
	}

	private static PSMData getPeptide(int index) {
		double precursorMZ;
		byte precursorCharge;
		String peptideModSeq;
		float retentionTime;
		if (index==0) {
			precursorMZ=500.730213;
			precursorCharge=(byte)2;
			peptideModSeq="MQS[+80.0]LSLNK";
			retentionTime=1198.3428f;
		} else if (index==1) {
			precursorMZ=500.899664;
			precursorCharge=(byte)3;
			peptideModSeq="SRPTS[+80.0]FADELAAR";
			retentionTime=1591.183f;
		} else if (index==2) {
			precursorMZ=503.551374;
			precursorCharge=(byte)3;
			peptideModSeq="A[+42.0]QRHS[+80.0]DSSLEEK";
			retentionTime=517.8737f;
		} else if (index==3) {
			precursorMZ=503.75254;
			precursorCharge=(byte)4;
			peptideModSeq="ASQEPS[+80.0]PKPGTEVIPAAPR";
			retentionTime=1328.8877f;
		} else if (index==4) {
			precursorMZ=767.052219;
			precursorCharge=(byte)3;
			peptideModSeq="DKRPLS[+80.0]GPDVGTPQPAGLASGAK";
			retentionTime=3686.5938f;
		} else if (index==10) {
			precursorMZ=492.241695;
			precursorCharge=(byte)3;
			peptideModSeq="KTAPTLS[+80.0]PEHWK";
			retentionTime=3873.0986f;
		} else {
			precursorMZ=503.272853;
			precursorCharge=(byte)3;
			peptideModSeq="KLS[+80.0]SGDLRVPVTR";
			retentionTime=1309.1414f;
		}
		
		PSMData psmdata=new PSMData(new HashSet<String>(), 0, precursorMZ, precursorCharge, peptideModSeq, retentionTime, 0, 0, 12);
		return psmdata;
	}
	
	public void testGetUniqueFragmentIonsRightInwards() {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		String peptide="LYSGS[+80.0]PTR";
		byte precursorCharge=2;
		
		ArrayList<String> peptideModSeqs=PhosphoPermuter.getPermutations(peptide, params.getAAConstants());

		HashMap<String, FragmentationModel> entryMap=new HashMap<String, FragmentationModel>();
		for (String peptideModSeq : peptideModSeqs) {
			FragmentationModel model=new FragmentationModel(peptideModSeq, params.getAAConstants());
			entryMap.put(peptideModSeq, model);
		}

		System.out.println("RIGHT TO LEFT");
		// go right to left, drop the first
		for (int i=peptideModSeqs.size()-1; i>=1; i--) {
			String targetPeptide=peptideModSeqs.get(i);
			AmbiguousPeptideModSeq targetPeptideName=AmbiguousPeptideModSeq.getRightAmbiguity(targetPeptide, AmbiguousPeptideModSeq.modifiableAAs, params.getAAConstants());
			
			HashMap<String, FragmentationModel> modelBatch=new HashMap<String, FragmentationModel>();
			// shrink the number of unique ions subtractors to the pool of remaining sequences to the right
			for (int j=0; j<=i; j++) {
				String seq=peptideModSeqs.get(j);
				modelBatch.put(seq, entryMap.get(seq));
			}
			System.out.println(targetPeptide+" SIZE: "+modelBatch.size());
			FragmentIon[] targets=PhosphoLocalizer.getUniqueFragmentIons(targetPeptide, precursorCharge, modelBatch, params);
			System.out.println(targetPeptideName.getPeptideAnnotation()+": "+General.toString(targets));
			System.out.println(targetPeptideName);
			System.out.println();
			assertTrue(targets.length<=16);
		}

		System.out.println("LEFT TO RIGHT");
		// go left to right, drop the last
		for (int i=0; i<peptideModSeqs.size()-1; i++) {
			String targetPeptide=peptideModSeqs.get(i);
			AmbiguousPeptideModSeq targetPeptideName=AmbiguousPeptideModSeq.getLeftAmbiguity(targetPeptide, AmbiguousPeptideModSeq.modifiableAAs, params.getAAConstants());

			HashMap<String, FragmentationModel> modelBatch=new HashMap<String, FragmentationModel>();
			// shrink the number of unique ions subtractors to the pool of
			// remaining sequences to the right
			for (int j=peptideModSeqs.size()-1; j>=i; j--) {
				String seq=peptideModSeqs.get(j);
				modelBatch.put(seq, entryMap.get(seq));
			}
			System.out.println(targetPeptide+" SIZE: "+modelBatch.size());
			FragmentIon[] targets=PhosphoLocalizer.getUniqueFragmentIons(targetPeptide, precursorCharge, modelBatch, params);
			System.out.println(targetPeptideName.getPeptideAnnotation()+": "+General.toString(targets));
			System.out.println(targetPeptideName);
			System.out.println();
			assertTrue(targets.length<=16);
		}
	}
	
	public void testGetUniqueFragmentIons() {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		
		byte precursorCharge=1;		
		String peptide="S[+80.0]SSSR";
		
		ArrayList<String> peptideModSeqs=PhosphoPermuter.getPermutations(peptide, params.getAAConstants());

		HashMap<String, FragmentationModel> entryMap=new HashMap<String, FragmentationModel>();
		for (String peptideModSeq : peptideModSeqs) {
			FragmentationModel model=new FragmentationModel(peptideModSeq, params.getAAConstants());
			entryMap.put(peptideModSeq, model);
		}
		
		for (int i=0; i<peptideModSeqs.size()-1; i++) {
			String targetPeptide=peptideModSeqs.get(i);
			
			HashMap<String, FragmentationModel> modelBatch=new HashMap<String, FragmentationModel>();
			// shrink the number of unique ions subtractors to the pool of remaining sequences to the right
			for (int j=peptideModSeqs.size()-1; j>=i; j--) {
				String seq=peptideModSeqs.get(j);
				modelBatch.put(seq, entryMap.get(seq));
			}
			FragmentIon[] uniqueIons=PhosphoLocalizer.getUniqueFragmentIons(targetPeptide, precursorCharge, modelBatch, params);
			assertTrue(uniqueIons.length<10);
		}
	}

	public void testIonGeneration() {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		
		byte charge=1;		
		String peptide="S[+80.0]SSR";
		
		ArrayList<String> permutations=PhosphoPermuter.getPermutations(peptide, params.getAAConstants());
		assertEquals(3, permutations.size());
		assertTrue(permutations.contains("S[+79.96633]SSR"));
		assertTrue(permutations.contains("SS[+79.96633]SR"));
		assertTrue(permutations.contains("SSS[+79.96633]R"));
		
		HashMap<String, FragmentationModel> entryMap=new HashMap<String, FragmentationModel>();
		for (String peptideModSeq : permutations) {
			FragmentationModel model=new FragmentationModel(peptideModSeq, params.getAAConstants());
			entryMap.put(peptideModSeq, model);
			FragmentIon[] ions=model.getPrimaryIonObjects(params.getFragType(), charge);
			//System.out.println(peptideModSeq);
			for (int i=0; i<ions.length; i++) {
				//System.out.println(ions[i]+"\t"+ions[i].mass);
			}
			//System.out.println();
		}

		HashMap<String, FragmentIon[]> uniqueIons=PhosphoLocalizer.getUniqueFragmentIons(charge, entryMap, params);
		
		assertEquals("b1,y3", General.toString(uniqueIons.get("S[+79.96633]SSR")));
		assertEquals("b2,y2", General.toString(uniqueIons.get("SSS[+79.96633]R")));
	}
}
