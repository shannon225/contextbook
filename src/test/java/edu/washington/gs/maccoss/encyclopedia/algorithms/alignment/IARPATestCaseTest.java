package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.LibraryReportExtractor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideReportData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.QuantitativeDIAData;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class IARPATestCaseTest {
	public static TObjectFloatHashMap<String> getGlassSlideTargetMap() {
		TObjectFloatHashMap<String> targetMap=new TObjectFloatHashMap<>();
		targetMap.put("GQHSSGSGQSPGHGQR",100000f);
		targetMap.put("FPSVSLQEASSFFR",150000f);
		targetMap.put("HSASQEGQDTIR",200000f);
		targetMap.put("QGSSAGSSSSYGQHGSGSR",200000f);
		targetMap.put("VPVDVAYR",200000f);
		targetMap.put("DILTIDISR",250000f);
		targetMap.put("DILTIDIGR",300000f);
		targetMap.put("VQYDLQK",300000f);
		targetMap.put("TYLISSIPLHGAFNYK",300000f);
		targetMap.put("FGQGVHHAAGQAGNEAGR",400000f);
		targetMap.put("AEGPEVDVNLPK",500000f);
		targetMap.put("AIGGALSSVGGGSSTIK",500000f);
		targetMap.put("HSGIGHGQASSAVR",500000f);
		targetMap.put("AEAEALYQTK",1000000f);
		targetMap.put("AGGSYGFGGAR",1000000f);
		targetMap.put("AIGGGLSSVGGGSSTIK",1000000f);
		targetMap.put("ELHPVLK",1000000f);
		targetMap.put("QNLEPLFEQYINNLR",1000000f);
		targetMap.put("SSGGSSSVK",1000000f);
		targetMap.put("SWNGSVEILK",1000000f);
		targetMap.put("TYLISSIPLQGAFNYK",1000000f);
		targetMap.put("YGQHGSGSR",1000000f);
		targetMap.put("YQELQITAGR",1000000f);
		targetMap.put("ALETVQER",2000000f);
		return targetMap;
	}
	
	public static TObjectFloatHashMap<String> getTargetMap() {
		TObjectFloatHashMap<String> targetMap=new TObjectFloatHashMap<>();
		targetMap.put("AEAEALYQIK", 5000000f);
		targetMap.put("AEAEALYQTK", 1000000f);
		targetMap.put("QVPGFGVADALGNR", 100000f);
		targetMap.put("AEGPEVDVNLPK", 150000f);
		targetMap.put("AGGSYGFGGAR", 200000f);
		targetMap.put("HAGIGHGQASSAVR", 200000f);
		targetMap.put("HSASQEGQDTIR", 200000f);
		targetMap.put("LVWEDTLDK", 250000f);
		targetMap.put("FPSVSLQEASSFFR", 300000f);
		targetMap.put("FGQGVHHAAAQAGNEAGR", 400000f);
		targetMap.put("YATTAYVPSEEINLVVK", 400000f);
		targetMap.put("GSGLGAGQGTNGASVK", 500000f);
		targetMap.put("SALSGHLETVILGLLK", 500000f);
		targetMap.put("SISVSVAGGALLGR", 500000f);
		targetMap.put("GETISGGNFHGEYPAK", 800000f);
		targetMap.put("AIGGGLSSVGGGSSTIK", 1000000f);
		targetMap.put("DILTIDIGR", 1000000f);
		targetMap.put("DVTVLQNTDGNNNEAWAK", 1000000f);
		targetMap.put("ELHPVLK", 1000000f);
		targetMap.put("FGQGVHHAAGQAGNEAGR", 1000000f);
		targetMap.put("GILIDTSR", 1000000f);
		targetMap.put("GQHSSGSGQSPGHDQR", 1000000f);
		targetMap.put("GQHSSGSGQSPGHGQR", 1000000f);
		targetMap.put("GSGLGAGQGSNGASVK", 1000000f);
		targetMap.put("HGLVATHTLTVR", 1000000f);
		targetMap.put("IMDVHDGK", 1000000f);
		targetMap.put("LVWEDTLVK", 1000000f);
		targetMap.put("QGSSAGSSSSYGQHGSGSR", 1000000f);
		targetMap.put("QNLEPLFEQYINNLR", 1000000f);
		targetMap.put("SSGGSSSVK", 1000000f);
		targetMap.put("VPVDVAYR", 1000000f);
		targetMap.put("YTQTFTLHANR", 1000000f);
		targetMap.put("TYLISSIPLQGAFNYK", 1200000f);
		targetMap.put("DILTIDISR", 1500000f);
		targetMap.put("HSGIGHGQASSAVR", 1500000f);
		targetMap.put("PEPSISLEPR", 1500000f);
		targetMap.put("GETVSGGNFHGEYPAK", 2000000f);
		targetMap.put("QCAELETAIADAEQR", 2000000f);
		targetMap.put("YGQHGSGSR", 2000000f);
		targetMap.put("NWNGSVEILK", 4000000f);
		targetMap.put("ALETVQER", 5000000f);
		targetMap.put("SWNGSVEILK", 5000000f);
		targetMap.put("YGQHGSGSCQSSGHGR", 5000000f);
		targetMap.put("AIGGALSSVGGGSSTIK", 8000000f);
		targetMap.put("QGSGSGQSPGHGQR", 10000000f);
		targetMap.put("VQYDLQK", 10000000f);
		targetMap.put("QGSSAGSSSSYGPHGSGSR", 20000000f);
		return targetMap;
	}

	
	private static final PecanSearchParameters PARAMETERS = new PecanSearchParameters(new AminoAcidConstants(),
			FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"),
			false, true, false);


	public static void main(String[] args) throws Exception {
		//File referenceFile = new File("/Volumes/searle_ssd/iarpa/final_stormy/individuals/final_stormy_individuals_clib.elib");
		//TObjectFloatHashMap<String> targetMap=getTargetMap();

		File referenceFile = new File("/Volumes/searle_ssd/iarpa/final_stormy/glassSlides/glass_slides_clib.elib");
		TObjectFloatHashMap<String> targetMap=getGlassSlideTargetMap();
		
		LibraryFile library=new LibraryFile();
		library.openFile(referenceFile);
		
		Pair<ArrayList<String>, ArrayList<PeptideReportData>> pair=LibraryReportExtractor.extractMatrix(library, PARAMETERS.getAAConstants());
		
		
		ArrayList<String> samples=pair.x;
		System.out.print("Peptide");
		for (String string : samples) {
			System.out.print("\t"+string);
		}
		System.out.println();
		for (PeptideReportData data : pair.y) {
			if (targetMap.contains(data.getPeptideSeq())) {
				System.out.print(data.getPeptideModSeq());
				
				float[] totalTICs=new float[samples.size()];
				
				for (int i=0; i<totalTICs.length; i++) {
					QuantitativeDIAData quantitativeData=data.getQuantitativeData(samples.get(i));

					if (quantitativeData!=null) {
						totalTICs[i]=quantitativeData.getTIC();
					}
				}
				
				for (int i=0; i<totalTICs.length; i++) {
					if (totalTICs[i]>targetMap.get(data.getPeptideSeq())) {
						System.out.print("\t"+totalTICs[i]);
					} else {
						System.out.print("\t"+0.0f);
					}
				}
				System.out.println();
			}
		}
	}
}
