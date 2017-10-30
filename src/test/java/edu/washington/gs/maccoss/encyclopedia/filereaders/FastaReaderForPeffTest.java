package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.ExtendedFastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TCharDoubleHashMap;
import junit.framework.TestCase;

public class FastaReaderForPeffTest extends TestCase {
	private final static int minLength = 8;
	private final static int maxLength = 40;
	private final static int maxMissedCleavages = 1;
	private final static AminoAcidConstants constants = new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());
	private final static double mzMin = 400d;
	private final static double mzMax = 900d;
	private final static double windowSize = 10d;
	
	
	public static void main(String[] args) throws Exception {
		File peffFile=new File("J:/1_LabData/20171017_peff_fileformat/nextprot2017_testPEFF1.0rc25_a.peff");
		//File peffFile=new File("J:/1_LabData/20171017_peff_fileformat/nextprot2017_testPEFF1.0rc25_small.peff");
		// File outputFile
		InputStream is=new FileInputStream(peffFile);
		long startTime=System.currentTimeMillis();
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(new BufferedReader(new InputStreamReader(is)), peffFile.getName(), null, true);
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
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
			/*
			for (String peptide : enzyme.digestProtein(entry, minLength, maxLength, maxMissedCleavages, constants)) {				
			}
			*/
			count++;
		}
		endTime=System.currentTimeMillis();
		System.out.println("\n"+entries.size()+" entries\t"+sum+" variants");
		duration=(endTime-startTime);
		System.out.println("total time: "+duration+" ms");
		
		

		is=new FileInputStream(peffFile);
		startTime=System.currentTimeMillis();
		entries=FastaReader.readFasta(new BufferedReader(new InputStreamReader(is)), peffFile.getName(), null, false);
		enzyme=DigestionEnzyme.getEnzyme("trypsin");
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
			/*
			for (String peptide : enzyme.digestProtein(entry, minLength, maxLength, maxMissedCleavages, constants)) {				
			}
			*/
			count++;
		}
		endTime=System.currentTimeMillis();
		System.out.println("\n"+entries.size()+" entries\t"+sum+" variants");
		duration=(endTime-startTime);
		System.out.println("total time: "+duration+" ms");

		
	}
	
	private static void compute(ArrayList<String> peptides, FastaEntryInterface entry) {
		
	}
	
	public void testFunction(){
		
	}
	
}

