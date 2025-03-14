package edu.washington.gs.maccoss.encyclopedia.gui.dia.interactive;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.Optional;
import java.util.StringTokenizer;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.HasRetentionTime;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideAccessionMatchingTrie;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SimplePeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;

public class InteractivePeptidePrecursor extends SimplePeptidePrecursor implements HasRetentionTime {
	private static final AminoAcidConstants aaConstants=new AminoAcidConstants();
	private final float rtInSecs;
	private final double precursorMZ;
	
	private byte isPassing; // CAN BE MUTATED!
	private Range rtRangeInSecs; // CAN BE MUTATED!
	
	public static void main(String[] args) throws Exception {
		File chromatogrindrCSV=new File("/Users/searleb/Downloads/Chimerys_Chromatogrindr_IL2.csv");
		File fastaFile=new File("/Users/searleb/Downloads/thorium/mus_musculus_reviewed_uniprot.fasta");
		File newLibraryFile=new File("/Users/searleb/Downloads/Chimerys_Chromatogrindr_IL2.dlib");
		
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
	
		final ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>(); 
		TableParser.parseCSV(chromatogrindrCSV, new TableParserMuscle() {
			
			@Override
			public void processRow(Map<String, String> row) {
				String peptideModSeq=row.get("Peptide");
				byte charge=Byte.parseByte(row.get("Charge"));
				float rtInMin=Float.parseFloat(row.get("RT (min)"));
				
				double mz=params.getAAConstants().getChargedMass(peptideModSeq, charge);
				int group=(int)Math.floor(mz/100.0);
				if (group<4) group=4;
				if (group>9) group=9;
				
				HashSet<String> accessions=new HashSet<String>();
				AnnotatedLibraryEntry entry=FragmentationModel.generateEntry(peptideModSeq, "GFP_"+group, accessions, charge, rtInMin*60f, false, params);
				entries.add(entry);
			}
			
			@Override
			public void cleanup() {
			}
		});

		ArrayList<FastaEntryInterface> proteins=FastaReader.readFasta(fastaFile, params);
		PeptideAccessionMatchingTrie trie=new PeptideAccessionMatchingTrie(entries, Optional.ofNullable(params.getEnzyme()));
		trie.addFasta(proteins);

		LibraryFile newLibrary=new LibraryFile();
		newLibrary.openFile();
		Logger.logLine("Writing library file "+newLibrary.getName());
		newLibrary.dropIndices();
		newLibrary.addEntries(entries);
		newLibrary.addProteinsFromEntries(entries);
		newLibrary.createIndices();
		newLibrary.saveAsFile(newLibraryFile);
	}
	
	public InteractivePeptidePrecursor(LibraryEntry entry) {
		this(entry.getPeptideModSeq(), entry.getPrecursorCharge(), entry.getRetentionTimeInSec());
	}

	public InteractivePeptidePrecursor(String peptideModSeq, byte precursorCharge, float rtInSecs) {
		this(peptideModSeq, precursorCharge, rtInSecs, (byte)0);
	}

	public InteractivePeptidePrecursor(String peptideModSeq, byte precursorCharge, float rtInSecs, byte isPassing) {
		this(peptideModSeq, precursorCharge, rtInSecs, isPassing, null);
	}

	public InteractivePeptidePrecursor(String peptideModSeq, byte precursorCharge, float rtInSecs, byte isPassing, Range rtRangeInSecs) {
		super(peptideModSeq, precursorCharge, aaConstants);
		this.rtInSecs=rtInSecs;
		this.precursorMZ=aaConstants.getChargedMass(peptideModSeq, precursorCharge);
		this.isPassing=isPassing;
		this.rtRangeInSecs=rtRangeInSecs;
	}
	
	public static String getHeader() {
		return "#Peptide\tRT(min)\tCharge\tIsPassing\tRTStart(min)\tRTStop(min)";
	}
	
	public String toString() {
		StringBuilder sb=new StringBuilder();
		sb.append(getPeptideModSeq());
		sb.append('\t');
		sb.append(rtInSecs/60f);
		sb.append('\t');
		sb.append(getPrecursorCharge());
		sb.append('\t');
		sb.append(getIsPassing());
		
		if (rtRangeInSecs==null) {
			sb.append('\t');
			sb.append('\t');
		} else {
			sb.append('\t');
			sb.append(rtRangeInSecs.getStart()/60f);
			sb.append('\t');
			sb.append(rtRangeInSecs.getStop()/60f);
		}
		return sb.toString();
	}
	
	public static InteractivePeptidePrecursor parseFromLine(String line) {
		StringTokenizer st = new StringTokenizer(line);
		String peptideModSeq = st.nextToken();
		float rtMin = Float.parseFloat(st.nextToken());
		byte charge = Byte.parseByte(st.nextToken());
		byte passes = st.hasMoreTokens() ? Byte.parseByte(st.nextToken()) : 0;
		
		Range range=null;
		if (st.hasMoreTokens()) {
			float first =Float.parseFloat(st.nextToken());
			float second =Float.parseFloat(st.nextToken());
			range=new Range(first*60f, second*60f);
		}
		peptideModSeq=PeptideUtils.getCorrectedMasses(peptideModSeq);
		return new InteractivePeptidePrecursor(peptideModSeq, charge, rtMin * 60f, passes, range);
	}

	@Override
	public float getRetentionTimeInSec() {
		return rtInSecs;
	}
	
	public double getPrecursorMZ() {
		return precursorMZ;
	}
	
	/**
	 * -1 for not passing
	 *  0 for unassigned
	 *  1 for passing
	 * @return
	 */
	public byte getIsPassing() {
		return isPassing;
	}

	/**
	 * -1 for not passing
	 *  0 for unassigned
	 *  1 for passing
	 */
	public void setIsPassing(boolean pass) {
		if (pass) {
			isPassing=1;
		} else {
			isPassing=-1;
		}
	}
	public void removeIsPassing() {
		isPassing=0;
	}
	
	/**
	 * CAN BE NULL
	 * @return
	 */
	public Range getRTRange() {
		return rtRangeInSecs;
	}
	public void setRtRangeInSecs(Range rtRangeInSecs) {
		this.rtRangeInSecs=rtRangeInSecs;
	}
}
