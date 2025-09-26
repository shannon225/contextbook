package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.util.CombinatoricsUtils;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.ExtendedFastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filewriters.FastaWriter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassErrorUnitType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.procedure.TIntProcedure;
import gnu.trove.map.hash.TCharIntHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.procedure.TCharIntProcedure;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

public class FastaReaderTest extends TestCase {


	public static void mainS(String[] args) throws Exception {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		File f=new File("/Users/searleb/Downloads/thorium/mus_musculus_reviewed_uniprot.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, params);
		String[] targets = new String[] { "P04202", "Q5SVQ0", "Q69ZK0", "P06240", "P43404", "A2AJK6", "P48453",
				"Q06180", "P29352", "P06800", "P27870", "Q61160", "P10417", "P03958", "P17095", "O08900", "O70494",
				"P06332", "Q08775", "P97305", "Q9JKD8", "P47930", "Q03267", "P08046", "Q8BZZ3", "Q61735", "O54824",
				"Q9Z1X4", "P10417", "P04187", "Q00417", "Q61164", "P04202", "P67778", "Q80TQ2", "P34902", "Q08024",
				"P63158", "Q8K1X4", "P01590", "Q99JB2", "P09055", "P70224", "Q60766", "Q8BWF2", "Q99MI6", "Q02111",
				"P58801", "Q3TTA7", "Q8K352", "O35129", "Q80V62", "Q99JG7", "Q8CIS0", "Q8R0W6", "Q60611", "Q64131",
				"O08573", "P42082", "Q60974", "Q03347", "P15379", "P22646", "Q8C863", "Q924T7", "Q3U2S4", "Q80TQ2",
				"Q08024", "Q8K1X4", "Q91XB0", "P70224", "Q9QZL0", "P29352", "Q61160", "P10417", "Q60611", "O54839",
				"Q64131", "P14429", "P06339", "P01897", "P01898", "P01901", "Q03347", "P10300", "P01731", "P15379",
				"P24788", "P37217", "Q78PG9", "Q8CGF7", "P16297", "P01590" };
		boolean[] found=new boolean[targets.length];

		ArrayList<FastaEntryInterface> saved=new ArrayList<FastaEntryInterface>();
		for (FastaEntryInterface entry : entries) {
			for (int i = 0; i < targets.length; i++) {
				if (entry.getAccession().indexOf(targets[i])>=0) {
					saved.add(entry);
					
					System.out.println(targets[i]+"\t"+params.getEnzyme().digestProtein(entry, 7, 40, 0, params.getAAConstants(), false).size());
					found[i]=true;
					break;
				}
			}
		}
		System.out.println(saved.size());
		for (int i = 0; i < found.length; i++) {
			if (!found[i]) System.out.println(targets[i]);
		}

		FastaWriter.writeFasta(new File("/Users/searleb/Downloads/thorium/selected_mus_musculus_reviewed_uniprot.fasta"), saved);
	}
	
	public static void mainW(String[] args) throws Exception {
		PecanSearchParameters parameters=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), false, true, false);
		File f=new File("/Users/searleb/Documents/maccoss/barnes/chris_barnes/mus_musculus_reviewed_uniprot.fasta"); 
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);

		for (FastaEntryInterface entry : entries) {
			ArrayList<FastaPeptideEntry> peptides=parameters.getEnzyme().digestProtein(entry, parameters.getMinPeptideLength(), parameters.getMaxPeptideLength(), parameters.getMaxMissedCleavages(), parameters.getAAConstants(), parameters.isRequireVariableMods());
			for (FastaPeptideEntry peptide : peptides) {
				for (byte z = 1; z <=4; z++) {
					double mass=parameters.getAAConstants().getChargedMass(peptide.getSequence(), z);
					if (mass>642.3&&mass<642.4) {
						char[] aas=peptide.getSequence().toCharArray();
						if (aas[2]=='N') {
						if (aas[3]=='I'||aas[3]=='L') {
								if (aas[4]=='I'||aas[4]=='L') {
									System.out.println(peptide.getSequence()+", "+z+", "+mass);
								}
							}
						}
					}
				}
			}
		}
	}

	public static void mainX(String[] args) throws Exception {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		File f=new File("/Users/searleb/Downloads/uniprotkb_proteome_UP000005640_2024_01_07.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, params);
		
		float[] ppms=new float[] {100000, 50000, 20000, 10000, 5000, 2000, 1000, 500, 200, 100, 50, 20, 10, 5, 2, 1, 0.5f, 0.2f, 0.1f, 0.05f, 0.02f, 0.01f};
		for (float ppm : ppms) {
			MassTolerance tolerance=new MassTolerance(ppm, MassErrorUnitType.PPM);
			
			TDoubleArrayList masses=new TDoubleArrayList();
			for (FastaEntryInterface entry : entries) {
				float mass=(float)params.getAAConstants().getMass(entry.getSequence());
				masses.add(mass);
			}
			
			masses.sort();
			int unique=0;
			for (int i = 0; i < masses.size(); i++) {
				double prev=i<=0?0.0:masses.get(i-1);
				double next=i>=masses.size()-1?Double.MAX_VALUE:masses.get(i+1);
				double curr=masses.get(i);
				
				if (tolerance.compareTo(next, curr)!=0&&tolerance.compareTo(prev, curr)!=0) {
					unique++;
				}
			}
			System.out.println(unique+" / "+entries.size()+" = "+(unique/(float)entries.size())+" at "+ppm+" ppm");
		}
	}

	public static void mainN(String[] args) throws Exception {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		File f=new File("/Users/searleb/Downloads/uniprotkb_proteome_UP000005640_2024_01_07.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, params);
		
		int[] density=new int[4000]; // 1000 Da per bin
		for (FastaEntryInterface entry : entries) {
			float mass=(float)params.getAAConstants().getMass(entry.getSequence());
			
			
			int index=Math.round((mass)/1000f);
			density[index]++;
		}
		
		for (int i = 1; i < density.length; i++) {
			System.out.println(Log.log10(i*1000)+"\t"+density[i]);
		}
	}
	
	// check for equal peptides between fastas
	public static void main(String[] args) throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		
		
		File f1=new File("/Users/searleb/Downloads/UP000001940_6239.fasta");
		File f2=new File("/Users/searleb/Downloads/Fasta_Target_only_CHIMERYS_Benchmark_human-canonical.fasta");
		DigestionEnzyme enzyme = DigestionEnzyme.getEnzyme("Trypsin");

		HashSet<String> matches=new HashSet<String>();
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f1, parameters);
		for (FastaEntryInterface entry : entries) {
			ArrayList<FastaPeptideEntry> peptides=enzyme.digestProtein(entry, 8, 40, 2, parameters.getAAConstants(), false);
			for (FastaPeptideEntry peptide : peptides) {
				matches.add(peptide.getSequence().replace("I", "L"));
			}			
		}

		HashSet<String> matches2=new HashSet<String>();
		ArrayList<FastaEntryInterface> entries2=FastaReader.readFasta(f2, parameters);
		for (FastaEntryInterface entry : entries2) {
			ArrayList<FastaPeptideEntry> peptides=enzyme.digestProtein(entry, 8, 40, 2, parameters.getAAConstants(), false);
			for (FastaPeptideEntry peptide : peptides) {
				matches2.add(peptide.getSequence().replace("I", "L"));
			}			
		}
		
		int collisions=0;
		for (String peptide : matches) {
			if (matches2.contains(peptide)) {
				collisions++;
			}
		}
		
		System.out.println(collisions+", "+matches.size()+", "+matches2.size());		
	}
	
	// check for reverse and forward matches
	public static void mainJ(String[] args) throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		File f=new File("/Users/searleb/Downloads/UP000001940_6239.fasta");
		DigestionEnzyme enzyme = DigestionEnzyme.getEnzyme("Trypsin");

		int counter=0;
		int total=0;
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);
		for (FastaEntryInterface entry : entries) {
			ArrayList<FastaPeptideEntry> peptides=enzyme.digestProtein(entry, 7, 40, 2, parameters.getAAConstants(), false);
			for (FastaPeptideEntry peptide : peptides) {
				total++;
				String inner=peptide.getSequence().substring(1, peptide.getSequence().length()-1);
				String reverse=new String(General.reverse(inner.toCharArray()));
				if (inner.equals(reverse)) {
					counter++;
				}
			}
			
		}
		
		System.out.println(counter+" / "+total+" of "+entries.size()+" proteins --> "+(100*counter/(float)total));
	}
	
	public static void mainGelAccuracy(String[] args) throws Exception {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		File f=new File("/Users/searleb/Downloads/uniprotkb_proteome_UP000005640_2024_01_07.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, params);
		float error=0.05f; //(+/-5%)
		
		int[] density=new int[40000]; // 100 Da per bin
		for (FastaEntryInterface entry : entries) {
			float mass=(float)params.getAAConstants().getMass(entry.getSequence());
			
			int minIndex=Math.round((mass-mass*error)/100f);
			int maxIndex=Math.round((mass+mass*error)/100f);
			for (int i = minIndex; i <= maxIndex; i++) {
				if (i>0&&i<density.length) {
					density[i]++;
				}
			}
		}
		
		for (int i = 1; i < density.length; i++) {
			System.out.println(Log.log10(i*100)+"\t"+density[i]);
		}
	}
	
	public static void mainL(String[] args) throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		String fastaFilePath = "/Users/searleb/Documents/OSU/teaching/proteomics_class/homework/homework3/UP000000625_83333.fasta";
		DigestionEnzyme enzyme = DigestionEnzyme.getEnzyme("Trypsin");
		
		File f=new File(fastaFilePath);

		int proteinCount=0;
		int previousDec=0;
		int counter=0;
		TCharIntHashMap aaGlobalMap=new TCharIntHashMap();
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);
		for (FastaEntryInterface entry : entries) {
			proteinCount++;
			int thisDec = Math.round(100f*proteinCount/(entries.size()));
			if (thisDec>previousDec) {
				previousDec=thisDec;
				//System.out.println(thisDec*1+"%");
			}
			ArrayList<FastaPeptideEntry> peptides=enzyme.digestProtein(entry, 7, 40, 2, parameters.getAAConstants(), false);
			counter+=peptides.size();

			TCharIntHashMap aaMap=new TCharIntHashMap();
			
			for (char c : entry.getSequence().toCharArray()) {
				aaMap.adjustOrPutValue(c, 1, 1);
				aaGlobalMap.adjustOrPutValue(c, 1, 1);
			}


//			System.out.print(entry.getAccession());
//			for (char c : AminoAcidConstants.AAs) {
//				System.out.print("\t"+aaMap.get(c));
//			}
//			System.out.println();
		}
		//System.out.println(counter);
		

		for (char c : AminoAcidConstants.AAs) {
			System.out.print("\t"+aaGlobalMap.get(c));
		}
	}
	
	public static void main10(String[] args) throws Exception {
		PecanSearchParameters parameters=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), false, true, false);
		
		LibraryFile library=new LibraryFile();
		library.openFile(new File ("/Users/searleb/Documents/teaching/encyclopedia/hela_pecan_narrow_library_604to616.dlib"));
		
		ArrayList<LibraryEntry> spectra=library.getAllEntries(false, parameters.getAAConstants());
		HashSet<String> accessions=new HashSet<>();
		for (LibraryEntry entry : spectra) {
			accessions.addAll(entry.getAccessions());
		}
		
		File f=new File("/Users/searleb/Documents/teaching/encyclopedia/uniprot-9606.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);
		FastaWriter writer=new FastaWriter(new File ("/Users/searleb/Documents/teaching/encyclopedia/uniprot-9606_604to616.fasta"));
		for (FastaEntryInterface entry : entries) {
			if (accessions.contains(entry.getAccession())) {
				writer.write(entry);
			}
		}
		writer.close();
	}
	
	public static void mainT(String[] args) throws Exception {
		PecanSearchParameters parameters=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), false, true, false);
		byte charge=6;
		File f=new File("/Users/searleb/Downloads/20201207_smallLibrary-1.fasta"); 
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);
		
		for (FastaEntryInterface entry : entries) {
			System.out.println(parameters.getAAConstants().getChargedMass(entry.getSequence(), charge));
		}
	}

	/**
	 * counts the number of accessions per peptide
	 * @param args
	 * @throws Exception
	 */
	public static void mainA(String[] args) throws Exception {
		PecanSearchParameters parameters=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), false, true, false);
		File f=new File("/Users/searleb/Downloads/bo_files/dmel-all-translation-r5.57_biognosysiRT_WR.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);
		
		Range[] ranges=new Range[] {new Range(350, 1250)};
		
		HashMap<String, HashSet<String>> accessions=new HashMap<>();
		
		for (Range range : ranges) {
			for (FastaEntryInterface entry : entries) {
				ArrayList<FastaPeptideEntry> peptides=parameters.getEnzyme().digestProtein(entry, parameters.getMinPeptideLength(), parameters.getMaxPeptideLength(), parameters.getMaxMissedCleavages(), parameters.getAAConstants(), parameters.isRequireVariableMods());
	
				for (FastaPeptideEntry peptide : peptides) {
					for (byte charge=parameters.getMinCharge(); charge<=parameters.getMaxCharge(); charge++) {
						double mz=parameters.getAAConstants().getChargedMass(peptide.getSequence(), charge);
	
						if (range.contains((float)mz)) {
							HashSet<String> pepAcc=accessions.get(peptide.getSequence());
							if (pepAcc==null) {
								pepAcc=new HashSet<>();
								accessions.put(peptide.getSequence(), pepAcc);
							}
							pepAcc.addAll(peptide.getAccessions());
						}
					}
				}
			}
		}

		int[] counts=new int[21];
		for (HashSet<String> entry : accessions.values()) {
			int size=Math.min(counts.length-1, entry.size());
			counts[size]++;
		}
		Logger.logLine("Accession count histogram: ");
		for (int i=0; i<counts.length; i++) {
			Logger.logLine(i+" Acc\t"+counts[i]+" Counts");
		}
	}
	
	public static void mainQ(String[] args) throws Exception {
		PecanSearchParameters parameters=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), false, true, false);
		File f=new File("/Users/searleb/Downloads/swissprot_reviewed_9606_13dec2019_20379_entries.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);
		
		Range[] ranges=new Range[] {new Range(400, 500), new Range(600, 700), new Range(900, 1000)};
		
		for (Range range : ranges) {
			int count=0;
			for (FastaEntryInterface entry : entries) {
				ArrayList<FastaPeptideEntry> peptides=parameters.getEnzyme().digestProtein(entry, parameters.getMinPeptideLength(), parameters.getMaxPeptideLength(), parameters.getMaxMissedCleavages(), parameters.getAAConstants(), parameters.isRequireVariableMods());
	
				for (FastaPeptideEntry peptide : peptides) {
					for (byte charge=parameters.getMinCharge(); charge<=parameters.getMaxCharge(); charge++) {
						double mz=parameters.getAAConstants().getChargedMass(peptide.getSequence(), charge);
	
						if (range.contains((float)mz)) {
							count++;
						}
					}
				}
			}
			System.out.println(range+" --> "+count);
		}
	}
	/*public static void main(String[] args) throws IOException, FileNotFoundException {
		File f=new File("/Users/searleb/Downloads/hg38_6FT.fasta");
		File out=new File("/Users/searleb/Downloads/hg38_coding.fasta");
		BufferedReader in=new BufferedReader(new FileReader(f));
		FastaWriter writer=new FastaWriter(out);

		boolean inAnnotation=false;
		StringBuilder annotation=new StringBuilder();
		StringBuilder sequence=new StringBuilder();
		char[] buffer=new char[1024*1024];
		int length;
		int counter=0;
		int index=0;
		int startIndex=0;
		while ((length=in.read(buffer))>=0) {
			counter++;
			if (counter%100==0) {
				System.out.println(" "+counter+" MB");
			} else if (counter%10==0) {
				System.out.print(". ");
			} else {
				System.out.print('.');
			}
			
			for (int i=0; i<length; i++) {
				index++;
				if (buffer[i]=='>') {
					inAnnotation=true;
					annotation.setLength(0);
				} else if (buffer[i]=='\n') {
					if (sequence.length()>=8) {
						writer.write(new FastaEntry(f.getName(), annotation.toString()+"."+(startIndex+1)+"."+(index-1), sequence.toString()));
					}
					sequence.setLength(0);
					
					inAnnotation=false;
					index=0;
					startIndex=index;
				} else {
					if (inAnnotation) {
						annotation.append(buffer[i]);
					} else {
						if (buffer[i]=='*'||buffer[i]=='X') {
							if (sequence.length()>=8) {
								writer.write(new FastaEntry(f.getName(), annotation.toString()+"."+(startIndex+1)+"."+(index-1), sequence.toString()));
							}
							sequence.setLength(0);
							startIndex=index;
						} else {
							sequence.append(buffer[i]);
						}
					}
				}
			}
		}
		in.close();
		writer.close();
	}*/
	
	
	/**
	 * for pecan peptide listing
	 * @param args
	 * @throws Exception
	 */
	public static void main2(String[] args) throws Exception {
		PecanSearchParameters parameters=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), false, true, false);
		//File f=new File("/Users/searleb/Documents/projects/phosphopedia/sp_iso_HUMAN_4.9.2015_UP000005640.fasta");
		File f=new File("/Users/searleb/Documents/chromatogram_library_manuscript/real_pecan/cerevisiae_orf_trans_all.fasta");
		ArrayList<FastaEntryInterface> targetProteins=FastaReader.readFasta(f, parameters);
		
		//PrintWriter writer=new PrintWriter("/Users/searleb/Documents/chromatogram_library_manuscript/sp_iso_HUMAN_4.9.2015_UP000005640.peptides.txt");
		PrintWriter writer=new PrintWriter("/Users/searleb/Documents/chromatogram_library_manuscript/real_pecan/cerevisiae_orf_trans_all.peptides.txt");
		for (FastaEntryInterface entry : targetProteins) {
			ArrayList<FastaPeptideEntry> peptides=parameters.getEnzyme().digestProtein(entry, parameters.getMinPeptideLength(), parameters.getMaxPeptideLength(), parameters.getMaxMissedCleavages(), parameters.getAAConstants(), false);
			for (FastaPeptideEntry peptide : peptides) {
				writer.println(entry.getAccession()+"\t"+PeptideUtils.getPeptideSeq(peptide.getSequence()));
			}
		}
		writer.flush();
		writer.close();
	}
	
	public static void mainx(String[] args) throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		File f=new File("/Users/bsearle/Documents/prosit/hela/uniprot-9606.fasta");
		
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);
		
		for (FastaEntryInterface entry : entries) {
			System.out.println(entry.getSequence().length());
		}
	}
	
	public static void main4(String[] args) {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		//File f=new File("/Users/searleb/Documents/projects/phosphopedia/sp_iso_HUMAN_4.9.2015_UP000005640.fasta");
		File f=new File("/Users/searleb/Documents/school/uniprot-9606.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);
		AminoAcidConstants constants=new AminoAcidConstants();
		System.out.println(entries.size());

		int countKR=0;
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		for (FastaEntryInterface entry : entries) {
			int charge=1+getCount(entry.getSequence(), 'K', 'R');
			double mass=constants.getMass(entry.getSequence())+MassConstants.oh2;
			double chargedMass=(mass+MassConstants.protonMass*charge)/charge;
			
			//System.out.println(charge);
			ArrayList<FastaPeptideEntry> peptides=enzyme.digestProtein(entry, 8, 40, 2, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
			for (FastaPeptideEntry string : peptides) {
				int pepCharge=2;
				double pepMass=constants.getMass(string.getSequence())+MassConstants.oh2;
				double pepChargedMass=(pepMass+MassConstants.protonMass*pepCharge)/pepCharge;
				
				if (pepChargedMass>(665.3204-5)&&pepChargedMass<(665.3204+5)) {
					System.out.println(string);
				}
			}
		}
	}

	public static void main3(String[] args) {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		AminoAcidConstants aaConstants=parameters.getAAConstants();
		//File f=new File("/Users/searleb/Documents/projects/phosphopedia/sp_iso_HUMAN_4.9.2015_UP000005640.fasta");
		File f=new File("/Users/searleb/Documents/school/uniprot-9606.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);
		Range possiblePrecursors=new Range(400, 1000);
		
		DigestionEnzyme trypsin=DigestionEnzyme.getEnzyme("trypsin");
		DigestionEnzyme gluc=DigestionEnzyme.getEnzyme("glu-c");
		DigestionEnzyme chymotrypsin=DigestionEnzyme.getEnzyme("chymotrypsin");
		DigestionEnzyme lysc=DigestionEnzyme.getEnzyme("lys-c");
		DigestionEnzyme argc=DigestionEnzyme.getEnzyme("arg-c");
		DigestionEnzyme aspn=DigestionEnzyme.getEnzyme("asp-n");
			
		int totalTotalSTY=0;
		int totalVisibleSTY=0;
		int totalTheOnlyOneSTY=0;
		int totalOneOfTwoSTY=0;
		for (FastaEntryInterface entry : entries) {
			int totalSTY=getCount(entry.getSequence(), 'S', 'T', 'Y');
			TIntHashSet visibleSTYIndicies=new TIntHashSet();
			TIntHashSet theOnlyOneSTYIndicies=new TIntHashSet();
			TIntHashSet oneOfTwoSTYIndicies=new TIntHashSet();

			for (DigestionEnzyme enzyme : new DigestionEnzyme[] {trypsin, chymotrypsin}) { //,trypsin,chymotrypsin,gluc,aspn,lysc,argc
				ArrayList<FastaPeptideEntry> peptides=enzyme.digestProtein(entry, 0, 9999999, 0, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
				int peptideIndex=0;
				for (FastaPeptideEntry peptide : peptides) {
					String sequence=peptide.getSequence();
	
					double mass=aaConstants.getMass(sequence)+MassConstants.oh2+79.966331;
					boolean ok=false;
					for (int charge=2; charge<=4; charge++) {
						double chargedMass=(mass+MassConstants.protonMass*charge)/charge;
						if (possiblePrecursors.contains(chargedMass)) {
							ok=true;
							break;
						}
					}
	
					int countSTY=getCount(sequence, 'S', 'T', 'Y');
					
					if (ok) {
						visibleSTYIndicies.addAll(General.add(getAllIndicies(sequence, 'S', 'T', 'Y'), peptideIndex));
						if (countSTY==1) {
							theOnlyOneSTYIndicies.add(peptideIndex+getFirstIndex(sequence, 'S', 'T', 'Y'));
						}
						if (countSTY<=2) {
							oneOfTwoSTYIndicies.addAll(General.add(getAllIndicies(sequence, 'S', 'T', 'Y'), peptideIndex));
						}
					}
					
					peptideIndex+=sequence.length();
				}
			}
			totalTotalSTY+=totalSTY;
			totalVisibleSTY+=visibleSTYIndicies.size();
			totalTheOnlyOneSTY+=theOnlyOneSTYIndicies.size();
			totalOneOfTwoSTY+=oneOfTwoSTYIndicies.size();
		}
		System.out.println(totalTotalSTY+"\t"+totalVisibleSTY+"\t"+totalTheOnlyOneSTY+"\t"+totalOneOfTwoSTY);
	}
	
	public static void main5(String[] args) {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		//File f=new File("/Users/searleb/Documents/projects/phosphopedia/sp_iso_HUMAN_4.9.2015_UP000005640.fasta");
		File f=new File("/Users/searleb/Documents/school/uniprot-9606.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);

		int countBase=0;
		int countNTermProtein=0;
		int countNTermPyroGlu=0;
		int countTryp=0;
		int countMet=0;
		int countSTY=0;
		int countQN=0;
		int countKR=0;
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		for (FastaEntryInterface entry : entries) {
			countNTermProtein++;
			ArrayList<FastaPeptideEntry> peptides=enzyme.digestProtein(entry, 8, 40, 2, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
			for (FastaPeptideEntry peptide : peptides) {
				String sequence=peptide.getSequence();
				countBase++;
				countNTermProtein++;
				if (sequence.charAt(0)=='Q'||sequence.charAt(0)=='C') {
					countNTermPyroGlu+=2;
				} else {
					countNTermPyroGlu++;
				}
				
				countTryp+=getCombinatorial(sequence, 'W');
				countMet+=getCombinatorial(sequence, 'M');
				countSTY+=getCombinatorial(sequence, 'S', 'T', 'Y');
				countQN+=getCombinatorial(sequence, 'Q', 'N');
				countKR+=getCombinatorial(sequence, 'K', 'R');
			}
		}
		System.out.println(countBase+"\tcountBase");
		System.out.println(countNTermProtein+"\tcountNTermProtein");
		System.out.println(countNTermPyroGlu+"\tcountNTermPyroGlu");
		System.out.println(countTryp+"\tcountTryp");
		System.out.println(countMet+"\tcountMet");
		System.out.println(countQN+"\tcountQN");
		System.out.println(countSTY+"\tcountSTY");
	}
	
	static int getCombinatorial(String sequence, char... target) {
		int num=getCount(sequence, target);
		if (num==0) return 1;
		if (num==1) return 2;
		if (num==2) return 4;
		return (int)(1+num+CombinatoricsUtils.binomialCoefficient(num, 2)+CombinatoricsUtils.binomialCoefficient(num, 3));
	}

	private static int getFirstIndex(String sequence, char... target) {
		int lowestIndex=Integer.MAX_VALUE;
		
		for (int i=0; i<target.length; i++) {
			int index=sequence.indexOf(target[i]);
			if (index>=0&&lowestIndex>index) {
				lowestIndex=index;
			}
		}
		return lowestIndex;
	}

	private static int[] getAllIndicies(String sequence, char... target) {
		TIntArrayList indicies=new TIntArrayList();
		for (int index=0; index<sequence.length(); index++) {
			for (int i=0; i<target.length; i++) {
				if (sequence.charAt(index)==target[i]) {
					indicies.add(index);
				}
			}
		}
		return indicies.toArray();
	}

	private static int getCount(String sequence, char... target) {
		int num=0;
		for (char c : sequence.toCharArray()) {
			for (int i=0; i<target.length; i++) {
				if (c==target[i]) {
					num++;
				}
			}
		}
		return num;
	}
	
	public void testFastaParsing() {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		String bsa=">ALBU_HUMAN Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2\n"+"MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF\n"
				+"EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEP\n"+"ERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLF\n"
				+"FAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAV\n"+"ARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLK\n"
				+"ECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYAR\n"+"RHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFE\n"
				+"QLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVV\n"+"LNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTL\n"
				+"SEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLV\n"+"AASQAALGL";

		FastaEntryInterface entry=FastaReader.readFasta(bsa, "", parameters).get(0);
		assertEquals("ALBU_HUMAN", entry.getAccession());
		assertEquals("MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL", entry.getSequence());

		String ecoli=">gi|16131183|ref|NP_417763.1| 50S ribosomal subunit protein L18 [Escherichia coli str. K-12 substr. MG1655]\n"
				+"MDKKSARIRRATRARRKLQELGATRLVVHRTPRHIYAQVIAPNGSEVLVAASTVEKAIAEQLKYTGNKDA\n"+"AAAVGKAVAERALEKGIKDVSFDRSGFQYHGRVQALADAAREAGLQF\n";

		entry=FastaReader.readFasta(ecoli, "", parameters).get(0);
		assertEquals("MDKKSARIRRATRARRKLQELGATRLVVHRTPRHIYAQVIAPNGSEVLVAASTVEKAIAEQLKYTGNKDAAAAVGKAVAERALEKGIKDVSFDRSGFQYHGRVQALADAAREAGLQF", entry.getSequence());
		ArrayList<FastaPeptideEntry> peptides=DigestionEnzyme.getEnzyme("trypsin").digestProtein(entry, 8, 40, 1, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
		boolean found=false;
		for (FastaPeptideEntry fastaPeptideEntry : peptides) {
			if ("GIKDVSFDR".equals(fastaPeptideEntry.getSequence())) {
				found=true;
				break;
			}
		}
		assertTrue(found);
	}

	public void testFastaReader() throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		InputStream is=getClass().getResourceAsStream("/ecoli-190209-contam_correctNL.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(is, "ecoli-190209-contam_correctNL.fasta", parameters);

		assertEquals(4178, entries.size());
		
		File temp=File.createTempFile("ecoli_", ".fasta"); // test write
		temp.deleteOnExit();
		FastaWriter.writeFasta(temp, entries);
		ArrayList<FastaEntryInterface> writtenEntries=FastaReader.readFasta(temp, parameters);

		assertEquals(4178, writtenEntries.size());
		for (int i=0; i<writtenEntries.size(); i++) {
			assertEquals(entries.get(i).getAccession(), writtenEntries.get(i).getAccession());
			assertEquals(entries.get(i).getSequence(), writtenEntries.get(i).getSequence());
		}
		
		entries=FastaReader.readFasta(temp, "NP_", parameters); // test filter
		assertEquals(3934, entries.size());

		File reverseConcatenated=PercolatorExecutor.getFastaPlusDecoyFile(temp, parameters);
		ArrayList<FastaEntryInterface> reverseConcatenatedEntries=FastaReader.readFasta(reverseConcatenated, parameters);
		assertEquals(4178*2, reverseConcatenatedEntries.size());
		reverseConcatenated.deleteOnExit();

		/*
		TIntObjectHashMap<TFloatArrayList> peptideDefects=new TIntObjectHashMap<TFloatArrayList>();
		
		for (FastaEntry entry : entries) {
			ArrayList<String> peptides=PARAMETERS.getEnzyme().digestProtein(entry.getSequence(), PARAMETERS.getMinPeptideLength(), PARAMETERS.getMaxPeptideLength(), PARAMETERS.getMaxMissedCleavages());
			for (String sequence : peptides) {
				FragmentationModel model=PeptideUtils.getPeptideModel(sequence);
				double[] ions=model.getPrimaryIons(PARAMETERS.getFragType());
				for (double d : ions) {
					d=(d+1.00727646681290)/2;
					int nominalMass=(int)d;
					float defect=(float)(d-nominalMass);
					TFloatArrayList list=peptideDefects.get(nominalMass);
					if (list==null) {
						list=new TFloatArrayList();
						peptideDefects.put(nominalMass, list);
					}
					list.add(defect);
				}
			}
		}
		
		int[] keys=peptideDefects.keys();
		Arrays.sort(keys);
		for (int nominal : keys) {
			float[] defects=peptideDefects.get(nominal).toArray();
			Arrays.sort(defects);
			System.out.println(nominal+"\t"+defects[(int)(defects.length*0.05f)]+"\t"+defects[(int)(defects.length*0.25f)]+"\t"+defects[(int)(defects.length*0.5f)]+"\t"+defects[(int)(defects.length*0.75f)]+"\t"+defects[(int)(defects.length*0.95f)]);
		}
		*/
	}
	
	/*
	 * @MoMo
	 * test reading Peff format and using ExtendedFastaEntry 
	 */
	public void testFastaReaderForPeff() throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		InputStream is=getClass().getResourceAsStream("/nextprot2017_testPEFF1.0rc25_small.peff");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(new BufferedReader(new InputStreamReader(is)), "nextprot2017_testPEFF1.0rc25_small.peff", null, true, parameters);
		assertEquals(25, entries.size());

		for (FastaEntryInterface entry : entries) {
			if (!(entry instanceof ExtendedFastaEntry)) {
				throw new Exception("Error occured when reading peff file, each entry should be an ExtendedFastaEntry object");
			}
		}
	}

	public void testReadNonexisting() throws Exception {
		final Path tmpDir = Files.createTempDirectory(getName());
		FileUtils.forceDeleteOnExit(tmpDir.toFile());

		try {
			final Path fasta = Files.createTempFile(tmpDir, getName(), ".fasta");
			Files.delete(fasta);

			final ArrayList<FastaEntryInterface> entries = FastaReader.readFasta(fasta.toFile(), SearchParameterParser.getDefaultParametersObject());
			assertNotNull(entries);
			assertEquals("Didn't get zero entries from nonexistent fasta!",
					0, entries.size()
			);
		} finally {
			FileUtils.deleteQuietly(tmpDir.toFile());
		}
	}
}
