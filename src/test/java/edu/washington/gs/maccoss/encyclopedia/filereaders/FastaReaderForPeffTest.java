package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

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
	private final static int maxMissedCleavages = 0;
	private final static AminoAcidConstants constants = new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());

	public static void main(String[] args) throws Exception {
		
		
		
		//System.out.println(abc.substring(0, idx));
		
		/*
		File peffFile=new File("J:/1_LabData/20171017_peff_fileformat/nextprot2017_testPEFF1.0rc25_a.peff");
		//File outputFile
		InputStream is= new FileInputStream(peffFile);
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(new BufferedReader(new InputStreamReader(is)), peffFile.getName(), null, true);
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		System.out.println("Number of entry in peff file: "+entries.size());
		
		
		
		for (FastaEntryInterface entry : entries) {
			ArrayList<String> peptides=enzyme.digestProtein(entry, minLength, maxLength, maxMissedCleavages, constants);
			//
		}
		*/
		
	}
	
	public void testFunction(){
		
	}



}
