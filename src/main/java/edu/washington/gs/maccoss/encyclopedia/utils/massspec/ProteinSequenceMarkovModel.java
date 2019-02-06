package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import edu.washington.gs.maccoss.encyclopedia.algorithms.SSRCalc;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filewriters.FastaWriter;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearDiscriminantAnalysis;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedianDouble;
import gnu.trove.list.array.TIntArrayList;

public class ProteinSequenceMarkovModel {
	public static void main(String[] args) throws Exception {
		
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		File f=new File("/Users/bsearle/Documents/prosit/hela/uniprot-9606.fasta");
		//File f=new File("/Users/bsearle/Documents/prosit/ecoli/uniprot_ecoli_25jan2019.fasta");
		//File f=new File("/Users/bsearle/Documents/prosit/a_thaliana/uniprot_a_thaliana_30jan2019.fasta");
		
		File ann=new File("/Users/bsearle/Documents/prosit/hela/human_subcellular_localization.txt");
		HashMap<String, String> annotationMap=new HashMap<>();
		TableParser.parseTSV(ann, new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				annotationMap.put(row.get("id"), row.get("cc"));
			}
			
			@Override
			public void cleanup() {
			}
		});
		
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);
		
		ArrayList<String> annotations=new ArrayList<>();
		ArrayList<double[]> rows=new ArrayList<>();
		for (FastaEntryInterface entry : entries) {
			String sequence=entry.getSequence();
			
			rows.add(calculateStatistics(sequence));
			
			String accession=entry.getAccession();
			int lastIndexOf=accession.lastIndexOf('|');
			if (lastIndexOf>=0) {
				accession=accession.substring(lastIndexOf+1);
			}
			String a = annotationMap.get(accession);
			if (a==null) a="unknown";
			annotations.add(a);
		}

		double[][] transpose=new double[rows.get(0).length][];
		for (int i = 0; i < transpose.length; i++) {
			transpose[i]=new double[rows.size()];
		}
		for (int i = 0; i < rows.size(); i++) {
			double[] row=rows.get(i);
			for (int j = 0; j < row.length; j++) {
				transpose[j][i]=row[j];
			}
		}

		double[] min=new double[transpose.length];
		double[] max=new double[transpose.length];
		for (int i = 0; i < transpose.length; i++) {
			min[i]=QuickMedianDouble.select(transpose[i], 0.025);
			max[i]=QuickMedianDouble.select(transpose[i], 0.975);
		}
		double[] range=General.subtract(max, min);
		for (double[] row : rows) {
			// normalize column to 0-1
			for (int i = 0; i < row.length; i++) {
				if (row[i]<min[i]) {
					row[i]=0.0;
				} else if (row[i]>max[i]) {
					row[i]=1.0;
				} else {
					row[i]=(row[i]-min[i])/range[i];
				}
			}
			double rowSum=General.sum(row);

			// then normalize row to sum to 1
			for (int i = 0; i < row.length; i++) {
				row[i]=row[i]/rowSum;
			}
		}

		ArrayList<double[]> positiveData=new ArrayList<>();
		ArrayList<double[]> negativeData=new ArrayList<>();
		for (int i = 0; i < rows.size(); i++) {
			//System.out.println((annotations.get(i).indexOf("brane")>0)+"\t"+General.toString(rows.get(i), "\t"));
			if (annotations.get(i).indexOf("brane")>0) {
				positiveData.add(rows.get(i));
			} else {
				negativeData.add(rows.get(i));
			}
		}
		
		LinearDiscriminantAnalysis lda=LinearDiscriminantAnalysis.buildModel(positiveData.toArray(new double[positiveData.size()][]), negativeData.toArray(new double[negativeData.size()][]));
		System.out.println(General.toString(lda.getCoefficients())+", "+lda.getConstant());


		// keep statistics
		int[] aaCount=new int[26];
		int[] lengthCount=new int[100];
		
		ProteinSequenceMarkovModel membraneModel=new ProteinSequenceMarkovModel();
		ProteinSequenceMarkovModel nonMembraneModel=new ProteinSequenceMarkovModel();
		
		for (int i = 0; i < rows.size(); i++) {
			//System.out.println((annotations.get(i).indexOf("brane")>0)+","+isMembrane(rows.get(i))+","+lda.getScore(General.toFloatArray(rows.get(i))));
			String sequence=entries.get(i).getSequence();
			for (char c : sequence.toCharArray()) {
				aaCount[getIndex(c)]++;
			}
			int index=Math.min(lengthCount.length-1, sequence.length()/10);
			lengthCount[index]++;
			
			if (isMembrane(rows.get(i))) {
				membraneModel.addProtein(sequence);
			} else {
				nonMembraneModel.addProtein(sequence);
			}
		}
		membraneModel.finalizeModel();
		nonMembraneModel.finalizeModel();

		int dbs=100;
		boolean write=true;
		
		int[][] randomAACounts=new int[26][];
		for (int i = 0; i < randomAACounts.length; i++) {
			randomAACounts[i]=new int[dbs];
		}
		int[][] randomLengthCounts=new int[100][];
		for (int i = 0; i < randomLengthCounts.length; i++) {
			randomLengthCounts[i]=new int[dbs];
		}
		
		for (int d = 0; d < dbs; d++) {
			FastaWriter writer=write?new FastaWriter(new File(f.getAbsolutePath()+"_"+d+".rand")):null;
			for (int i = 0; i < membraneModel.getNumberOfProteins(); i++) {
				String random=membraneModel.generateRandomSequence();
				for (char c : random.toCharArray()) {
					randomAACounts[getIndex(c)][d]++;
				}
				int index=Math.min(randomLengthCounts.length-1, random.length()/10);
				randomLengthCounts[index][d]++;
				
				FastaEntry entry=new FastaEntry(null, "RANDOM_MEMBRANE_"+i, random);
				if (write) writer.write(entry);
			}
	
			for (int i = 0; i < nonMembraneModel.getNumberOfProteins(); i++) {
				String random=nonMembraneModel.generateRandomSequence();
				for (char c : random.toCharArray()) {
					randomAACounts[getIndex(c)][d]++;
				}
				int index=Math.min(randomLengthCounts.length-1, random.length()/10);
				randomLengthCounts[index][d]++;
				
				FastaEntry entry=new FastaEntry(null, "RANDOM_NONMEMBRANE_"+i, random);
				if (write) writer.write(entry);
			}
			if (write) writer.close();
		}
		
		for (int i = 0; i < aaCount.length; i++) {
			float[] floatArray = General.toFloatArray(randomAACounts[i]);
			System.out.println(getChar(i)+","+aaCount[i]+", "+QuickMedian.select(floatArray, 0.25f)+", "+QuickMedian.select(floatArray, 0.5f)+", "+QuickMedian.select(floatArray, 0.75f));
		}
//		for (int i = 0; i < lengthCount.length; i++) {
//			float[] floatArray = General.toFloatArray(randomLengthCounts[i]);
//			System.out.println((i*10)+", "+lengthCount[i]+", "+QuickMedian.select(floatArray, 0.25f)+", "+QuickMedian.select(floatArray, 0.5f)+", "+QuickMedian.select(floatArray, 0.75f));
//		}
		
		/*BarnesHutTSne tsne;
		boolean parallel = false;
		if (parallel) {
			tsne = new ParallelBHTsne();
		} else {
			tsne = new BHTSne();
		}
		TSneConfiguration config = TSneUtils.buildConfig(rows.toArray(new double[rows.size()][]), 2, 3, 20.0f, 2000);
		double[][] Y = tsne.tsne(config);
		
		ArrayList<XYPoint> membrane=new ArrayList<>();
		ArrayList<XYPoint> other=new ArrayList<>();
		for (int i = 0; i < Y.length; i++) {
			if (annotations.get(i).indexOf("brane")>0) {
				membrane.add(new XYPoint(Y[i][0], Y[i][1]));
			} else {
				other.add(new XYPoint(Y[i][0], Y[i][1]));
			}
			//System.out.println(+"\t"+General.toString(Y[i], "\t"));
		}
		Charter.launchChart("xAxis", "yAxis", true, new XYTrace(membrane, GraphType.tinypoint, "membrane"), new XYTrace(other, GraphType.tinypoint, "other"));
		*/
		
		//System.out.println(model.lengthModel.lengths.size()+" / "+model.nodes.length);
		
		//System.out.println(random);
		
		//writer.println("modified_sequence,collision_energy,precursor_charge");
		//writer.close();
	}
	
	public static boolean isMembrane(double[] scores) {
		return (scores[0]*17.9375+scores[1]*-0.875+scores[2]*6.0-6.855)>=0;
	}
	
	private static final AminoAcidConstants aaConstants=new AminoAcidConstants();
	public static double[] calculateStatistics(String sequence) {
		double hydrophobicity=SSRCalc.getHydrophobicity(sequence);
		double mass=0.0;
		double charge=0.0;
		
		for (char aa : sequence.toCharArray()) {
			mass+=aaConstants.getMass(aa);
			
			if ('R'==aa||'K'==aa) {
				charge++;
			} else if ('D'==aa||'E'==aa) {
				charge--;
			}
		}
		return new double[] {hydrophobicity, Log.log10(mass), charge};
	}
	
	private final LengthModel lengthModel=new LengthModel();
	private final HashMap<String, Node> nodeMap=new HashMap<>();
	private Node[] nodes;
	private int[] centerCounts;
	private int totalCenterCounts;
	private int[] startCounts;
	private int totalStartCounts;
	private int[] stopCounts;
	private int totalStopCounts;
	
	public ProteinSequenceMarkovModel() {
		
	}
	
	public int getNumberOfProteins() {
		return lengthModel.lengths.size();
	}
	
	private static int hashSize=3;
	public void addProtein(String sequence) {
		if (sequence.length()<=hashSize) return;
		
		int startIndex = hashSize+1;
		loop: for (int i = startIndex; i < sequence.length(); i++) {
			// this is an index of the original string, so watch out
			String key=sequence.substring(i-hashSize-1, i-1);
			if (isBadSeq(key)) {
				continue loop;
			}
			char c=sequence.charAt(i);
			if (isBadAA(c)) {
				continue loop;
			}
			
			increment(key, c);
			if (i==startIndex) {
				incrementStart(key);
			}
		}
		String key=sequence.substring(sequence.length()-hashSize-1, sequence.length()-1);
		if (!isBadSeq(key)) {
			incrementStop(key);
		}

		lengthModel.increment(sequence.length());
	}

	private boolean isBadSeq(String seq) {
		for (int j = 0; j < seq.length(); j++) {
			char charAt = seq.charAt(j);
			if (isBadAA(charAt)) {
				return true;
			}
		}
		return false;
	}
	private boolean isBadAA(char charAt) {
		if (!Character.isUpperCase(charAt)) {
			return true;
		}
		if ('B'==charAt) {
			return true;
		} else if ('J'==charAt) {
			return true;
		} else if ('O'==charAt) {
			return true;
		} else if ('U'==charAt) {
			return true;
		} else if ('X'==charAt) {
			return true;
		} else if ('Z'==charAt) {
			return true;
		}
		return false;
	}

	private void increment(String key, char c) {
		Node n=nodeMap.get(key);
		if (n==null) {
			n=new Node(key.charAt(0), key.charAt(1), key.charAt(2));
			// this conserves memory by letting Java garbage collect the full sequence string 
			nodeMap.put(new String(key.toCharArray()), n); 
		}
		n.increment(c);
	}
	private void incrementStart(String key) {
		Node n=nodeMap.get(key);
		if (n==null) {
			n=new Node(key.charAt(0), key.charAt(1), key.charAt(2));
			// this conserves memory by letting Java garbage collect the full sequence string 
			nodeMap.put(new String(key.toCharArray()), n); 
		}
		n.incrementStart();
	}
	private void incrementStop(String key) {
		Node n=nodeMap.get(key);
		if (n==null) {
			n=new Node(key.charAt(0), key.charAt(1), key.charAt(2));
			// this conserves memory by letting Java garbage collect the full sequence string 
			nodeMap.put(new String(key.toCharArray()), n); 
		}
		n.incrementStop();
	}
	
	private void finalizeModel() {
		nodes=nodeMap.values().toArray(new Node[nodeMap.size()]);
		Arrays.sort(nodes);
		centerCounts=new int[nodes.length];
		for (int i = 0; i < centerCounts.length; i++) {
			centerCounts[i]=nodes[i].totalCount;
		}
		totalCenterCounts=General.sum(centerCounts);

		startCounts=new int[nodes.length];
		for (int i = 0; i < startCounts.length; i++) {
			startCounts[i]=nodes[i].starts;
		}
		totalStartCounts=General.sum(startCounts);

		stopCounts=new int[nodes.length];
		for (int i = 0; i < stopCounts.length; i++) {
			stopCounts[i]=nodes[i].stops;
		}
		totalStopCounts=General.sum(stopCounts);
	}
	
	private Node getRandomCenterNode(double prob) {
		int sum=(int)(prob*totalCenterCounts);
		
		int runningTotal=0;
		for (int i = 0; i < centerCounts.length; i++) {
			runningTotal+=centerCounts[i];
			if (runningTotal>sum) {
				return nodes[i];
			}
		}
		throw new RuntimeException("Over the counting index! Check for bugs.");
	}
	
	private Node getRandomStartNode(double prob) {
		int sum=(int)(prob*totalStartCounts);
		
		int runningTotal=0;
		for (int i = 0; i < startCounts.length; i++) {
			runningTotal+=startCounts[i];
			if (runningTotal>sum) {
				return nodes[i];
			}
		}
		throw new RuntimeException("Over the counting index! Check for bugs.");
	}
	
	private Node getRandomStopNode(double prob) {
		int sum=(int)(prob*totalStopCounts);
		
		int runningTotal=0;
		for (int i = 0; i < stopCounts.length; i++) {
			runningTotal+=stopCounts[i];
			if (runningTotal>sum) {
				return nodes[i];
			}
		}
		throw new RuntimeException("Over the counting index! Check for bugs.");
	}
	
	public String generateRandomSequence() {
		int length=lengthModel.getLength(Math.random());
		
		Node n=getRandomStartNode(Math.random());
		StringBuilder sb=new StringBuilder(n.toString());
		
		int middleLength=length-hashSize-hashSize;
		for (int i = 0; i < middleLength; i++) {
			sb.append(getNextAminoAcid(sb.substring(sb.length()-hashSize)));
		}
		
		n=getRandomStopNode(Math.random());
		sb.append(n.toString());
		return sb.toString();
	}
	
	private char getNextAminoAcid(String key) {
		Node n=nodeMap.get(key);
		if (n==null) {
			n=getRandomCenterNode(Math.random());
		}
		return n.getRandomAA(Math.random());
	}
	
	private class Node implements Comparable<Node> {
		
		private final char a;
		private final char b;
		private final char c;
		private final int[] counts=new int[26];
		private int totalCount=0;
		private int starts=0;
		private int stops=0;
		private int hash;
		
		public Node(char a, char b, char c) {
			this.a = a;
			this.b = b;
			this.c = c;
		}
		
		public String toString() {
			return new String(new char[] {a, b, c});
		}
		
		@Override
		public int compareTo(Node o) {
			if (o==null) return 1;
			int r=Character.compare(a, o.a);
			if (r!=0) return r;
			r=Character.compare(b, o.b);
			if (r!=0) return r;
			r=Character.compare(c, o.c);
			if (r!=0) return r;
			return 0;
		}
		
		@Override
		public boolean equals(Object obj) {
			if (obj instanceof Node) {
				return compareTo((Node)obj)==0;
			} else {
				return false;
			}
		}
		
		@Override
		public int hashCode() {
			int h = hash;
			if (h == 0) {
				char val[] = new char[] { a, b, c };

				for (int i = 0; i < 3; i++) {
					h = 31 * h + val[i];
				}
				hash = h;
			}
			return h;
		}
		
		private void increment(char c) {
			totalCount++;
			counts[getIndex(c)]++;
		}
		
		private void incrementStart() {
			starts++;
		}
		
		private void incrementStop() {
			stops++;
		}
		
		private char getRandomAA(double prob) {
			int sum=(int)(prob*totalCount);
			
			int runningTotal=0;
			for (int i = 0; i < counts.length; i++) {
				runningTotal+=counts[i];
				if (runningTotal>sum) {
					return getChar(i);
				}
			}
			throw new RuntimeException("Over the counting index! Check for bugs.");
		}
	}
	
	private class LengthModel {
		private final TIntArrayList lengths=new TIntArrayList();
		
		public LengthModel() {
		}
		
		private void increment(int length) {
			lengths.add(length);
		}
		
		public int getLength(double prob) {
			return lengths.get((int)(lengths.size()*prob));
		}
	}
	
	private static int getIndex(char c) {
		assert(c>=65&&c<=90);
		return c-65;
	}
	
	private static char getChar(int index) {
		return (char)(index+65);
	}
}
