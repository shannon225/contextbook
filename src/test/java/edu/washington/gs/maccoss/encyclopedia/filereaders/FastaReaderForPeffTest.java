package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import org.apache.commons.lang3.text.WordUtils;

import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.AlleleVariant;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.ExtendedFastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filewriters.FastaWriter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;


public class FastaReaderForPeffTest {
	private final static int minLength=8;
	private final static int maxLength=40;
	private final static int maxMissedCleavages=0;
	private final static AminoAcidConstants constants=new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());
	private final static DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");

	public static void main(String[] args) throws Exception {
		//checkDigestionRunningTimeForPeff();
		
		
		File peffFile=new File("J:/1_LabData/20171017_peff_fileformat/nextprot2017_testPEFF1.0rc25_a.peff");
		//File peffFile=new File("J:/1_LabData/20171017_peff_fileformat/nextprot2017_testPEFF1.0rc25_small.peff");
		File rangeFile=new File("J:/1_LabData/20171017_peff_fileformat/MZranges.csv");
		String outputFolder = "J:/1_LabData/20171017_peff_fileformat/";
		
		InputStream is=new FileInputStream(peffFile);
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(new BufferedReader(new InputStreamReader(is)), peffFile.getName(), null, true);
		ArrayList<Range> ranges=readRangeFile(rangeFile);
		Collections.sort(ranges);

		ArrayList<CountingResult> results=new ArrayList<CountingResult>();
		for (int i=0; i<3; i++) {
			results.add(new CountingResult(ranges.size()));
		}

		
		for (int i=0; i<entries.size(); i++) {
			FastaEntryInterface entry=entries.get(i);
			ArrayList<String> peptideSequences=enzyme.digestProtein(entry.getSequence(), minLength, maxLength, maxMissedCleavages, constants, new ArrayList<AlleleVariant>());
			ArrayList<AlleleVariant> variants=new ArrayList<AlleleVariant>();

			if (entry instanceof ExtendedFastaEntry) {
				variants=((ExtendedFastaEntry)entry).getPotentialVariant();
			}

			
			ArrayList<Peptide> peptides=new ArrayList<Peptide>();
			for (String sequence : peptideSequences) {
				int start=entry.getSequence().indexOf(sequence)+1;
				int end=start+sequence.length()-1;
				Peptide peptide=new Peptide(sequence, PeptideUtils.getExpectedChargeState(sequence));
				peptide.setSequenceStartAndEnd(start, end);
				peptides.add(peptide);
			}

			// only variants locate at the peptide sequence are associated with peptide
			Collections.sort(variants);
			associateVariantWithPeptide(variants, peptides);
			
			for (int chargeOffset=-1; chargeOffset<=1; chargeOffset++) {
				int index=chargeOffset+1;
				countPeptideAndVariantInRanges(peptides, ranges, chargeOffset, results.get(index));
			}
		}
		outputReports(results,ranges,outputFolder);
		
	}
	
	
	private static ArrayList<Range> readRangeFile(File rangeFile) {
		BufferedReader in=null;
		ArrayList<Range> ranges=new ArrayList<Range>();
		try {
			in=new BufferedReader(new InputStreamReader(new FileInputStream(rangeFile)));
		} catch (FileNotFoundException fne) {
			fne.printStackTrace();
			System.out.println("Can not find range csv file");
		}

		try {
			String eachline;
			while ((eachline=in.readLine())!=null) {
				if ((eachline.trim().length()==0)||!Character.isDigit(eachline.charAt(0))) {
					continue;
				}
				String[] mzs=eachline.split(",");
				try {
					ranges.add(new Range(Float.parseFloat(mzs[0]), Float.parseFloat(mzs[1])));
				} catch (NumberFormatException ne) {
					System.out.println("Error on parsing mz range from range csv file ["+eachline+"]");
				}
			}
		} catch (IOException ioe) {
			System.out.println("I/O Error found reading range csv file");
		} finally {
			if (in!=null) {
				try {
					in.close();
				} catch (IOException ioe) {
					ioe.printStackTrace();
				}
			}
		}
		return ranges;
	}
	
	private static void associateVariantWithPeptide(ArrayList<AlleleVariant> variants, ArrayList<Peptide> peptides) {
		int varIndex=0;
		for (int i=0; i<peptides.size(); i++) {
			Peptide peptide=peptides.get(i);
			for (int j=varIndex; j<variants.size(); j++) {
				AlleleVariant variant=variants.get(j);
				if (variant.getStopSite()<peptide.getStartSite()) {
					varIndex=j;
					continue;
				} else if (variant.getStartSite()>peptide.getEndSite()) {
					break;
				}
				peptide.addVariant(variants.get(j));
			}
		}
	}
	
	private static void countPeptideAndVariantInRanges(ArrayList<Peptide> peptides, ArrayList<Range> ranges, int chargeOffset, CountingResult result) {
		for (Peptide peptide : peptides) {

			if (peptide.getVariants().size()==0) {
				result.peptideWithoutVariant+=1;
				continue;
			}

			byte charge=(byte)(peptide.getCharge()+chargeOffset);
			if (charge==(byte)0) {
				result.peptideWithoutCharge+=1;
				continue;
			}

			double mz=constants.getChargedMass(peptide.getSequence(), charge);
			if ((mz>=ranges.get(ranges.size()-1).getStop())||(mz<ranges.get(0).getStart())) {
				result.peptideOutsider+=1;
				continue;
			}

			int rangeIndex=0;
			for (int i=0; i<ranges.size(); i++) {
				Range range=ranges.get(i);
				if ((mz>=range.getStart())&&(mz<range.getStop())) {
					rangeIndex=i;
					break;
				}
			}
			result.peptideInRanges[rangeIndex]+=1;

			boolean hasVariantInSameRange=false;
			for (int i=0; i<peptide.getVariants().size(); i++) {
				AlleleVariant variant=peptide.getVariants().get(i);

				// get variant sequence
				String variantSequence="";
				int offset=variant.getStartSite()-peptide.getStartSite();
				if (offset>0) {
					variantSequence+=peptide.getSequence().substring(0, offset);
				}
				if (!variant.getNewSequence().equals("*")) {
					variantSequence+=variant.getNewSequence();
					offset=peptide.getEndSite()-variant.getStopSite();
					if (offset>0) {
						variantSequence+=peptide.getSequence().substring(variant.getStopSite()-peptide.getStartSite()+1);
					}
				}

				// assuming variant sequence has same charge state
				double varMZ=constants.getChargedMass(variantSequence, charge);
				double delta=varMZ-mz;
				result.mzDifferences.add(delta);

				// counting variant type
				String sequenceVariation=variant.getOriginalSequence()+"->"+variant.getNewSequence();
				if (!result.variantTypeCount.containsKey(sequenceVariation)) {
					VariantType vType;
					if (variant.getNewSequence().equals("*")) {
						vType=new VariantType(sequenceVariation, Double.MAX_VALUE);
					} else {
						double deltaMass = constants.getMass(variant.getNewSequence())-constants.getMass(variant.getOriginalSequence());
						vType=new VariantType(sequenceVariation, deltaMass);
					}
					result.variantTypeCount.put(sequenceVariation, vType);
				}
				result.variantTypeCount.get(sequenceVariation).count+=1;

				// assign mz region
				if ((varMZ>=ranges.get(rangeIndex).getStart())&&(varMZ<ranges.get(rangeIndex).getStop())) {
					result.variantInSameRange[rangeIndex]+=1;
					hasVariantInSameRange=true;
				} else {
					result.variantInDiffRange[rangeIndex]+=1;
				}
			}
			if (hasVariantInSameRange) {
				result.peptideHasVariantInSamRange+=1;
			}
		}
	}
	
	private static void outputReports(ArrayList<CountingResult> results,ArrayList<Range> ranges,String outputFolder){
		
		PrintWriter writer; 
		File outputFile=new File(outputFolder+"sequence_variant_region_counting_summary.csv");
		try {
			writer=new PrintWriter(outputFile, "UTF-8");
		} catch (FileNotFoundException|UnsupportedEncodingException e) {
			throw new EncyclopediaException("Error writing peff file peptide-variant counting summry!", e);
		}

		for (int i=0; i<results.size(); i++) {
			writer.println("expected charge stat offset:,"+(i-1));
			writer.println("number of peptide without variant:,"+results.get(i).peptideWithoutVariant);
			writer.println("number of peptide with no charge:,"+results.get(i).peptideWithoutCharge);
			writer.println("number of peptide outside of range:,"+results.get(i).peptideOutsider);
			writer.println("number of peptide in mz regions:,"+General.sum(results.get(i).peptideInRanges));
			writer.println("number of peptide has one or more variants in same mz region:,"+results.get(i).peptideHasVariantInSamRange);
			writer.println("number of peptide-variant pair in same mz region:,"+General.sum(results.get(i).variantInSameRange));
			writer.println("number of peptide-variant pair in different mz region:,"+General.sum(results.get(i).variantInDiffRange));
			writer.println();
			
			String filepath = outputFolder+"sequence_variant_delta_mass_histogram_for_charge_offset_"+(i-1)+".csv";
			outputHistogramFile(results.get(i).mzDifferences,0.5f,filepath);
			filepath = outputFolder+"sequence_variant_mz_region_distribution_for_charge_offset_"+(i-1)+".csv";
			outputVariantMZRegionDistribution(ranges,results.get(i),filepath);
			filepath = outputFolder+"sequence_variation_record_for_charge_offset_"+(i-1)+".csv";
			outputSequenceVarianceFrequency(results.get(i),filepath);
		}
		writer.close();
		System.out.println("done");
	}
	
	private static void outputHistogramFile(ArrayList<Double> deltaMZs ,double binSize,String outputFilepath){
		Collections.sort(deltaMZs);
		double minValue= Math.floor(deltaMZs.get(0));
		double maxValue= Math.ceil(deltaMZs.get(deltaMZs.size()-1));
		int arraylength = (int)Math.ceil((maxValue-minValue)/binSize)+1;
		double[] x = new double[arraylength];
		int[] y = new int[arraylength];
		
		int index=0;
		x[0] = minValue;
		for (int i =0;i<x.length-1;i++){
			x[i+1]=x[i]+binSize;
			for (int j=index;j<deltaMZs.size();j++){
				double mz = deltaMZs.get(j);
				if (mz>=x[i+1]){
					index=j;
					break;
				} else {
					y[i]+=1;
				}
			}
		}
		PrintWriter writer; 
		File outputFile = new File(outputFilepath);
		try{
			writer=new PrintWriter(outputFile, "UTF-8");
		} catch (FileNotFoundException|UnsupportedEncodingException e) {
			throw new EncyclopediaException("Error writing peff file peptide-variant mz delta result!", e);
		}
		writer.println("deltaMS bin,counts");
		for (int i =0;i<x.length;i++){
			writer.println(x[i]+","+y[i]);
		}
		writer.close();
	}
	
	private static void outputVariantMZRegionDistribution(ArrayList<Range> ranges, CountingResult result, String outputFilepath) {
		PrintWriter writer;
		File outputFile=new File(outputFilepath);
		try {
			writer=new PrintWriter(outputFile, "UTF-8");
		} catch (FileNotFoundException|UnsupportedEncodingException e) {
			throw new EncyclopediaException("Error writing peff file peptide-variant mz delta result!", e);
		}

		writer.println("peptide m/z region,number of peptide,number of variant in same region, number of variant in different region");

		DecimalFormat mzFormat=new DecimalFormat("##.00");
		for (int i=0; i<ranges.size(); i++) {
			writer.println(mzFormat.format(ranges.get(i).getStart())+"-"+mzFormat.format(ranges.get(i).getStop())+","+result.peptideInRanges[i]+","+result.variantInSameRange[i]+","
					+result.variantInDiffRange[i]);
		}
		writer.close();
	}
	
	private static void outputSequenceVarianceFrequency (CountingResult result, String outputFilepath){
		PrintWriter writer;
		File outputFile=new File(outputFilepath);
		try {
			writer=new PrintWriter(outputFile, "UTF-8");
		} catch (FileNotFoundException|UnsupportedEncodingException e) {
			throw new EncyclopediaException("Error writing peff file peptide-variant mz delta result!", e);
		}

		writer.println("sequence variation,delta mass,count");

		DecimalFormat mzFormat=new DecimalFormat("##.00");
		for (String sequenceVariance : result.variantTypeCount.keySet()){
			writer.println(sequenceVariance+","+result.variantTypeCount.get(sequenceVariance).deltaMass()+","+result.variantTypeCount.get(sequenceVariance).count);
		}
		writer.close();		
	}
	




	private static void checkDigestionRunningTimeForPeff() throws FileNotFoundException {
		File peffFile=new File("J:/1_LabData/20171017_peff_fileformat/nextprot2017_testPEFF1.0rc25_a.peff");
		// File peffFile=new
		// File("J:/1_LabData/20171017_peff_fileformat/nextprot2017_testPEFF1.0rc25_small.peff");
		// File outputFile
		InputStream is=new FileInputStream(peffFile);
		long startTime=System.currentTimeMillis();
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(new BufferedReader(new InputStreamReader(is)), peffFile.getName(), null, true);
		
		System.out.println("Number of entry in peff file: "+entries.size());
		long endTime=System.currentTimeMillis();
		long duration=(endTime-startTime);
		System.out.println("reading time: "+duration+" ms\taverage: "+duration/entries.size());
		int count=0;
		int sum=0;
		startTime=System.currentTimeMillis();
		for (FastaEntryInterface entry : entries) {
			if (entry instanceof ExtendedFastaEntry) {
				sum+=((ExtendedFastaEntry)entry).getPotentialVariant().size();
			}
			if (count%2000==0) {
				System.out.print(count+" ");
			}
			enzyme.digestProtein(entry, minLength, maxLength, maxMissedCleavages, constants);
			count++;
		}
		endTime=System.currentTimeMillis();
		System.out.println("\n"+entries.size()+" entries\t"+sum+" variants");
		duration=(endTime-startTime);
		System.out.println("total time: "+duration+" ms");

		is=new FileInputStream(peffFile);
		startTime=System.currentTimeMillis();
		entries=FastaReader.readFasta(new BufferedReader(new InputStreamReader(is)), peffFile.getName(), null, false);
		
		System.out.println("Number of entry in peff file: "+entries.size());
		endTime=System.currentTimeMillis();
		duration=(endTime-startTime);
		System.out.println("reading time: "+duration+" ms\taverage: "+duration/entries.size());
		count=0;
		sum=0;
		startTime=System.currentTimeMillis();
		for (FastaEntryInterface entry : entries) {
			if (entry instanceof ExtendedFastaEntry) {
				sum+=((ExtendedFastaEntry)entry).getPotentialVariant().size();
			}
			if (count%2000==0) {
				System.out.print(count+" ");
			}
			enzyme.digestProtein(entry, minLength, maxLength, maxMissedCleavages, constants);
			count++;
		}
		endTime=System.currentTimeMillis();
		System.out.println("\n"+entries.size()+" entries\t"+sum+" variants");
		duration=(endTime-startTime);
		System.out.println("total time: "+duration+" ms");
	}

}

class Peptide implements Comparable<Peptide> {
	private String sequence="";
	private int start=0;
	private int end=0;
	private byte expectedCharge=(byte)0;
	private final ArrayList<AlleleVariant> variants=new ArrayList<AlleleVariant>();
	

	public Peptide(String sequence, byte charge) {
		this.sequence=sequence;
		this.expectedCharge=charge;
	}

	public void addVariant(AlleleVariant variant) {
		variants.add(variant);
	}

	public ArrayList<AlleleVariant> getVariants() {
		return variants;
	}

	public String getSequence() {
		return sequence;
	}

	public byte getCharge() {
		return expectedCharge;
	}

	public void setSequenceStartAndEnd(int start, int end) {
		this.start=start;
		this.end=end;
	}

	public int getStartSite() {
		return start;
	}

	public int getEndSite() {
		return end;
	}
	@Override
	public int compareTo(Peptide o) {
		if (this.getStartSite()>o.getStartSite()) {
			return 1;
		}
		return 0;
	}

}

class CountingResult {
	public int peptideWithoutVariant=0;
	public int peptideOutsider=0;
	public int peptideHasVariantInSamRange=0;
	public int peptideWithoutCharge=0;
	public final int[] peptideInRanges;
	public final int[] variantInSameRange;
	public final int[] variantInDiffRange;
	public ArrayList<Double> mzDifferences=new ArrayList<Double>();
	public HashMap<String, VariantType> variantTypeCount=new HashMap<String, VariantType>();

	public CountingResult(int noOfMZRegions) {
		peptideInRanges=new int[noOfMZRegions];
		variantInSameRange=new int[noOfMZRegions];
		variantInDiffRange=new int[noOfMZRegions];
	}
	 

}

class VariantType {
	private final String sequenceVariant;
	private final double deltaMass;
	public int count;
	public VariantType(String sequenceVariant, double delta){
		this.sequenceVariant = sequenceVariant;
		this.deltaMass=delta;
	}
	public String getSequenceVariant(){
		return sequenceVariant;
	}
	public double deltaMass (){
		return deltaMass;
	}
}

