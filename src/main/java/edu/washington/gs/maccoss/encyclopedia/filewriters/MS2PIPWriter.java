package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.StringTokenizer;
import java.util.zip.DataFormatException;

import org.apache.commons.lang3.StringUtils;

import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SimplePeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class MS2PIPWriter {
	private static final String NO_PTMS_CODE = "-";
	private static final String PTM_DELIMINATOR = "|";
	private static final String DEAMIDATED = "Deamidated";
	private static final String ACETYL = "Acetyl";
	private static final String PHOSPHO_Y = "PhosphoY";
	private static final String PHOSPHO_T = "PhosphoT";
	private static final String PHOSPHO_S = "PhosphoS";
	private static final String CARBAMIDOMETHYL = "Carbamidomethyl";
	private static final String OXIDATION = "Oxidation";
	private static final AminoAcidConstants emptyAAConstants=new AminoAcidConstants();

	public static void main(String[] args) throws Exception {
		File libraryFile = new File("/Volumes/bcsbluessd/swaney/CID_vs_HCD_enzymes/cid.dlib");
		LibraryFile library=new LibraryFile();
		library.openFile(libraryFile);
		
		writeMS2PIP("/Volumes/bcsbluessd/swaney/CID_vs_HCD_enzymes/cid.peprec", library, false);
	}

	public static void writeMS2PIP(String csvFileName, File fasta, DigestionEnzyme enzyme, byte minCharge, byte maxCharge, int maxMissedCleavages, Range mzRange, boolean addDecoys) throws FileNotFoundException {
		String fileName = checkPEPRECName(csvFileName, fasta, enzyme);
		Logger.logLine("Starting to build MS2PIP PEPREC: "+fileName);

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(fasta, parameters);
		
		HashSet<PeptidePrecursor> allPeptides=new HashSet<>();
		
		for (FastaEntryInterface entry : entries) {
			ArrayList<FastaPeptideEntry> peptidesInProtein=enzyme.digestProtein(entry, 6, 100, maxMissedCleavages, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
			for (FastaPeptideEntry pep : peptidesInProtein) {
				for (byte pepCharge = minCharge; pepCharge <=maxCharge; pepCharge++) {
					String seq=pep.getSequence();
					double pepMass=emptyAAConstants.getMass(seq)+MassConstants.oh2;
					double pepChargedMass=(pepMass+MassConstants.protonMass*pepCharge)/pepCharge;

					if (mzRange.contains(pepChargedMass)) {
						if (seq.indexOf('B')>=0||seq.indexOf('J')>=0||seq.indexOf('O')>=0||seq.indexOf('U')>=0||seq.indexOf('X')>=0||seq.indexOf('Z')>=0||seq.indexOf('*')>=0) {
							continue;
						} else {
							allPeptides.add(new SimplePeptidePrecursor(seq, pepCharge, emptyAAConstants));
						}
					}
				}
			}
		}

		int total = writePEPRECFile(fileName, addDecoys, allPeptides);
		Logger.logLine("Finished writing "+total+" peptides to MS2PIP PEPREC!");
	}
	
	public static void writeMS2PIP(String csvFileName, LibraryFile library, boolean addDecoys) throws FileNotFoundException {
		String fileName = checkPEPRECName(csvFileName, library.getFile(), null);
		Logger.logLine("Starting to build MS2PIP PEPREC: "+fileName);
		AminoAcidConstants aminoAcidConstants = new AminoAcidConstants();

		HashSet<PeptidePrecursor> allPeptides=new HashSet<>();
		try {
			ArrayList<LibraryEntry> allEntries = library.getAllEntries(false, aminoAcidConstants);
			allPeptides.addAll(allEntries);
		} catch (IOException e) {
			throw new EncyclopediaException("Error parsing library", e);
		} catch (SQLException e) {
			throw new EncyclopediaException("Error parsing library", e);
		} catch (DataFormatException e) {
			throw new EncyclopediaException("Error parsing library", e);
		}
		
		int total = writePEPRECFile(fileName, addDecoys, allPeptides);
		Logger.logLine("Finished writing "+total+" peptides to MS2PIP PEPREC!");
	}

	static String checkPEPRECName(String csvFileName, File fasta, DigestionEnzyme enzyme) {
		String fileName;
		if (null==csvFileName||StringUtils.isBlank(csvFileName)) {
			String enzymeText=enzyme==null?"":("."+enzyme.getPercolatorName());
			fileName = fasta.getAbsolutePath() + enzymeText + ".peprec";
		} else {
			fileName = csvFileName;
		}
		return fileName;
	}
	
	private static int writePEPRECFile(String fileName, boolean addDecoys, HashSet<PeptidePrecursor> allPeptides) throws FileNotFoundException {
		PrintWriter writer=new PrintWriter(fileName);
		int total=writePeptides(allPeptides, false, writer);
		if (addDecoys) {
			total+=writePeptides(allPeptides, true, writer);
		}
		writer.close();
		return total;
	}

	private static int writePeptides(HashSet<PeptidePrecursor> allPeptides, boolean writeDecoysInstead, PrintWriter writer) {
		HashSet<String> writtenPeptides=new HashSet<>();
		HashSet<String> ptmCodes=new HashSet<>();
		
		writer.println("spec_id modifications peptide charge");
		int total=1;
		EACHPEPTIDE: for (PeptidePrecursor peptidePrecursor : allPeptides) {
			if (writeDecoysInstead) {
				peptidePrecursor=new SimplePeptidePrecursor(PeptideUtils.reverse(peptidePrecursor.getPeptideModSeq(), emptyAAConstants), peptidePrecursor.getPrecursorCharge(), emptyAAConstants);
			}
			
			// remove peptides that don't match MS2PIP limitations:
			String seq = peptidePrecursor.getPeptideSeq();
			if (seq.indexOf('B')>=0||seq.indexOf('J')>=0||seq.indexOf('O')>=0||seq.indexOf('U')>=0||seq.indexOf('X')>=0||seq.indexOf('Z')>=0||seq.indexOf('*')>=0) {
				continue;
			}
			
			StringBuilder sb=new StringBuilder("p"+total+" ");
			
			FragmentationModel model=PeptideUtils.getPeptideModel(peptidePrecursor.getPeptideModSeq(), emptyAAConstants);
			
			String[] aas=model.getAas();
			if(seq.length()!=aas.length) {
				Logger.errorLine("Error parsing peptide sequence "+peptidePrecursor.getPeptideModSeq());
				continue EACHPEPTIDE;
			}
			
			double[] modificationMasses=model.getModificationMasses();
			boolean first=true;
			for (int i = 0; i < aas.length; i++) {
				String tag;
				if (!first) {
					tag=PTM_DELIMINATOR+(i+1)+PTM_DELIMINATOR;
				} else {
					tag=(i+1)+PTM_DELIMINATOR;
				}
				
				char aa=seq.charAt(i);
				int mass=(int)Math.round(modificationMasses[i]);
				if (aa=='M'&&mass==16) {
					sb.append(tag+OXIDATION);
					ptmCodes.add(OXIDATION);
					first=false;
				} else if (aa=='C') { //&&mass==57) { // ASSUMES CAM!
					sb.append(tag+CARBAMIDOMETHYL);
					ptmCodes.add(CARBAMIDOMETHYL);
					first=false;
				} else if (aa=='S'&&mass==80) {
					sb.append(tag+PHOSPHO_S);
					ptmCodes.add(PHOSPHO_S);
					first=false;
				} else if (aa=='T'&&mass==80) {
					sb.append(tag+PHOSPHO_T);
					ptmCodes.add(PHOSPHO_T);
					first=false;
				} else if (aa=='Y'&&mass==80) {
					sb.append(tag+PHOSPHO_Y);
					ptmCodes.add(PHOSPHO_Y);
					first=false;
				} else if (mass==42) {
					if (i==0) {
						// n-term
						tag=(i)+PTM_DELIMINATOR;
					}
					sb.append(tag+ACETYL);
					ptmCodes.add(ACETYL);
					first=false;
				} else if (aa=='N'&&mass==1) {
					sb.append(tag+DEAMIDATED);
					ptmCodes.add(DEAMIDATED);
					first=false;
				} else if (mass!=0) {
					Logger.errorLine("Error parsing peptide modfication "+peptidePrecursor.getPeptideModSeq()+" (only Oxidation, Carbamidomethyl, Phospho[STY], Acetyl, and Deamidated allowed)");
					continue EACHPEPTIDE;
				}
			}
			if (first) {
				sb.append(NO_PTMS_CODE);
			}

			sb.append(' ');
			sb.append(seq);
			sb.append(' ');
			sb.append(peptidePrecursor.getPrecursorCharge());
			
			String line=sb.toString();
			if (writtenPeptides.contains(line)) {
				continue;
			} else {
				writtenPeptides.add(line);
			}
			writer.println(line);
			total++;
		}
		
		StringBuilder sb=new StringBuilder();
		if (ptmCodes.size()>0) {
			sb.append("PTM codes:\n");
			for (String code : ptmCodes) {
				sb.append(code+","+getPTMMass(code)+","+getPTMTarget(code)+"\n");
			}	
		}
		Logger.logLine(sb.toString());
		return total;
	}
	
	public static String getPeptideModSeq(String peptideSequence, String ptmCodes) throws EncyclopediaException {
		if (ptmCodes==null||ptmCodes.length()==0||NO_PTMS_CODE.equals(ptmCodes)) return peptideSequence;
		
		StringTokenizer st=new StringTokenizer(ptmCodes, PTM_DELIMINATOR);
		char[] aaChars=peptideSequence.toCharArray();
		String[] aas=new String[aaChars.length];
		for (int i = 0; i < aas.length; i++) {
			aas[i]=Character.toString(aaChars[i]);
		}
		
		while (st.hasMoreTokens()) {
			int index=Integer.parseInt(st.nextToken());
			if (index>0) {
				index=index-1;
			}
			String code=st.nextToken();
			double mass = getPTMMass(code);
			
			if (mass!=0.0) {
				aas[index]=aas[index]+"["+mass+"]";
			}
		}
		
		StringBuilder sb=new StringBuilder();
		for (int i = 0; i < aas.length; i++) {
			sb.append(aas[i]);
		}
		return sb.toString();
	}

	private static double getPTMMass(String code) {
		double mass=0.0;
		if (OXIDATION.equalsIgnoreCase(code)) {
			mass=15.994915;
		} else if (CARBAMIDOMETHYL.equalsIgnoreCase(code)) {
			mass=57.0214635;
		} else if (PHOSPHO_S.equalsIgnoreCase(code)) {
			mass=PeptideModification.phosphorylation.getMass();
		} else if (PHOSPHO_T.equalsIgnoreCase(code)) {
			mass=PeptideModification.phosphorylation.getMass();
		} else if (PHOSPHO_Y.equalsIgnoreCase(code)) {
			mass=PeptideModification.phosphorylation.getMass();
		} else if (ACETYL.equalsIgnoreCase(code)) {
			mass=42.010565;
		} else if (DEAMIDATED.equalsIgnoreCase(code)) {
			mass=MassConstants.oh2-MassConstants.nh3;
		} else {
			Logger.errorLine("Error parsing peptide modfication "+code+" (only Oxidation, Carbamidomethyl, Phospho[STY], Acetyl, and Deamidated allowed)");
			throw new EncyclopediaException("Unexpected PTM code: ["+code+"]");
		}
		return mass;
	}

	private static String getPTMTarget(String code) {
		if (OXIDATION.equalsIgnoreCase(code)) {
			return "M";
		} else if (CARBAMIDOMETHYL.equalsIgnoreCase(code)) {
			return "C";
		} else if (PHOSPHO_S.equalsIgnoreCase(code)) {
			return "S";
		} else if (PHOSPHO_T.equalsIgnoreCase(code)) {
			return "T";
		} else if (PHOSPHO_Y.equalsIgnoreCase(code)) {
			return "Y";
		} else if (ACETYL.equalsIgnoreCase(code)) {
			return "N-term";
		} else if (DEAMIDATED.equalsIgnoreCase(code)) {
			return "N";
		} else {
			Logger.errorLine("Error parsing peptide modfication "+code+" (only Oxidation, Carbamidomethyl, Phospho[STY], Acetyl, and Deamidated allowed)");
			throw new EncyclopediaException("Unexpected PTM code: ["+code+"]");
		}
	}
}
