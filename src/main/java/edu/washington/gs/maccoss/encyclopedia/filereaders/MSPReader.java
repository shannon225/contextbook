package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.StringTokenizer;

import edu.washington.gs.maccoss.encyclopedia.algorithms.SSRCalc;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideAccessionMatchingTrie;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.map.hash.TIntDoubleHashMap;

public class MSPReader {
	public static void main(String[] args) throws Exception {
		File mspFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/cptac2_human_hcd_selected.msp");
		File fastaFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/UP000005640_9606.fasta");
		File libraryFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/cptac2_human_hcd_selected.elib");
		
		convertMSP(mspFile, fastaFile, libraryFile, SearchParameterParser.getDefaultParametersObject());
	}
	public static void convertMSP(File mspFile, File fastaFile, SearchParameters parameters) throws IOException, SQLException, IllegalArgumentException {
		String absolutePath=mspFile.getAbsolutePath();
		File libraryFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.DLIB);
		convertMSP(mspFile, fastaFile, libraryFile, parameters);
	}

	public static void convertMSP(File mspFile, File fastaFile, File libraryFile, SearchParameters parameters) throws IOException, SQLException, IllegalArgumentException {
		Logger.logLine("Reading MSP file "+mspFile.getName());
		ArrayList<LibraryEntry> entries=readMSP(mspFile, false);
		Logger.logLine("Read "+entries.size()+" total entries");
		
		Logger.logLine("Reading Fasta file "+fastaFile.getName());
		ArrayList<FastaEntryInterface> proteins=FastaReader.readFasta(fastaFile, parameters);
		Logger.logLine("Read "+proteins.size()+" total proteins");

		Logger.logLine("Constructing trie from library peptides");
		PeptideAccessionMatchingTrie trie=new PeptideAccessionMatchingTrie(entries);
		trie.addFasta(proteins);
		
		int[] counts=new int[21];
		for (LibraryEntry entry : entries) {
			int size=Math.min(counts.length-1, entry.getAccessions().size());
			counts[size]++;
		}
		Logger.logLine("Accession count histogram: ");
		for (int i=0; i<counts.length; i++) {
			Logger.logLine(i+" Acc\t"+counts[i]+" Counts");
		}

		LibraryFile library=new LibraryFile();
		library.openFile();
		library.dropIndices();

		if (counts[0]>0) {
			Logger.errorLine(counts[0]+" library entries can't be linked to proteins! These entries will be dropped.");
		}
		Logger.logLine("Writing library file "+libraryFile.getName());
		int batchSize=Math.max(1, entries.size()/10);
		int start=0;
		int stop=batchSize;
		while (true) {
			if (start>=entries.size()) break;
			ArrayList<LibraryEntry> sublist=new ArrayList<LibraryEntry>(entries.subList(start, stop));
			library.addEntries(sublist);
			library.addProteinsFromEntries(sublist);
			start=stop;
			stop=Math.min(entries.size(), stop+batchSize);
			Logger.logLine((start*100/entries.size())+"%");
		}
		library.createIndices();
		library.saveAsFile(libraryFile);
		library.close();
	}
	
	public static ArrayList<LibraryEntry> readMSP(File f, boolean keepAccessions) throws IOException, IllegalArgumentException{
		BufferedReader in=null;
		try {
			in=new BufferedReader(new FileReader(f));
			return readMSP(in, f.getName(), keepAccessions);

		} finally {
			if (in!=null) {
				try {
					in.close();
				} catch (IOException ioe) {
					throw(ioe);
				}
			}
		}
	}
	
	public static ArrayList<LibraryEntry> readMSP(String s, String fileName, boolean keepAccessions) throws IOException, IllegalArgumentException {
		return readMSP(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(s.getBytes(StandardCharsets.UTF_8)))), fileName, keepAccessions);
	}
	
	public static ArrayList<LibraryEntry> readMSP(InputStream s, String fileName, boolean keepAccessions) throws IOException, IllegalArgumentException {
		return readMSP(new BufferedReader(new InputStreamReader(s)), fileName, keepAccessions);
	}
	
	public static ArrayList<LibraryEntry> readMSP(BufferedReader in, String fileName, boolean keepAccessions) throws IOException, IllegalArgumentException {
		//TODO: take in parameters
		final AminoAcidConstants aaConstants = new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());

		ArrayList<LibraryEntry> entryList=new ArrayList<LibraryEntry>();
		String eachline=null;
		try {
			
			String peptideModSeq=null;
			String fullname=null;
			String accession=null;
			String altName=null; //Issue 90: A fallback when the full name is not present in the file
			byte precursorCharge=0;
			double precursorMZ=0.0;
			float retentionTime=0.0f;
			float score=0.0f;
			ArrayList<Peak> peaks=new ArrayList<Peak>();
			OUTERLOOP: while ((eachline=in.readLine())!=null) {
				eachline=eachline.trim();
				if (eachline.length()==0) {
					continue OUTERLOOP;
				}
				
				if (eachline.startsWith("Name: ")) {
					altName = eachline.substring(6);
					if (peaks.size()>0) {
						Pair<double[], float[]> peakArrays=Peak.toArrays(peaks);
						HashSet<String> accessions=new HashSet<String>();
						if (accession!=null) accessions.add(accession);
						
						if (retentionTime==0.0f) Logger.errorLine("MSP parsing error: expected to find retention time for "+peptideModSeq+" but it was missing."); 
						if (precursorMZ==0.0) Logger.errorLine("MSP parsing error: expected to find precursor M/Z for "+peptideModSeq+" but it was missing."); 
						if (precursorCharge==(byte)0) Logger.errorLine("MSP parsing error: expected to find charge for "+peptideModSeq+" but it was missing."); 
						LibraryEntry entry=new LibraryEntry(fileName, accessions, precursorMZ, precursorCharge, peptideModSeq, 1, retentionTime, score, peakArrays.x, peakArrays.y, aaConstants);
						retentionTime=0.0f;
						precursorMZ=0.0;
						precursorCharge=0;
						
						entryList.add(entry);
						if (entryList.size()%10000==0) {
							Logger.logLine("Read "+entryList.size()+" entries...");
						}
						peaks.clear();
						fullname=null;
					}
				} else if (eachline.startsWith("Num peaks: ")||eachline.startsWith("NumPeaks: ")) {
					INNERLOOP: while ((eachline=in.readLine())!=null) {
						eachline=eachline.trim();
						if (eachline.length()==0) {
							break INNERLOOP;
						}
						if (eachline.startsWith("Name: ")) {
							// indicates that we didn't get a normal linebreak to finish the entry. Assume that we've started the next entry.
							altName = eachline.substring(6);
							if (peaks.size()>0) {
								Pair<double[], float[]> peakArrays=Peak.toArrays(peaks);
								HashSet<String> accessions=new HashSet<String>();
								if (accession!=null) accessions.add(accession);
								
								if (retentionTime==0.0f) Logger.errorLine("MSP parsing error: expected to find retention time for "+peptideModSeq+" but it was missing."); 
								if (precursorMZ==0.0) Logger.errorLine("MSP parsing error: expected to find precursor M/Z for "+peptideModSeq+" but it was missing."); 
								if (precursorCharge==(byte)0) Logger.errorLine("MSP parsing error: expected to find charge for "+peptideModSeq+" but it was missing."); 
								LibraryEntry entry=new LibraryEntry(fileName, accessions, precursorMZ, precursorCharge, peptideModSeq, 1, retentionTime, score, peakArrays.x, peakArrays.y, aaConstants);
								retentionTime=0.0f;
								precursorMZ=0.0;
								precursorCharge=0;
								
								entryList.add(entry);
								if (entryList.size()%10000==0) {
									Logger.logLine("Read "+entryList.size()+" entries...");
								}
								peaks.clear();
								fullname=null;
							}
							break INNERLOOP;
						}
						
						StringTokenizer st=new StringTokenizer(eachline);
						double mass=Double.parseDouble(st.nextToken());
						float intensity=Float.parseFloat(st.nextToken());
						peaks.add(new Peak(mass, intensity));
					}
				} else if (eachline.startsWith("FullName: ")) {
					fullname=eachline.substring(10);
				} else if (eachline.startsWith("Charge: ")) {
					precursorCharge=Byte.parseByte(eachline.substring(8));
				} else if (eachline.startsWith("PrecursorMZ: ")) {
					precursorMZ=Double.parseDouble(eachline.substring(13));
				} else if (eachline.startsWith("RetentionTimeMins: ")) {
					retentionTime=Float.parseFloat(eachline.substring(19))*60f;
				} else if (eachline.startsWith("Comment: ")) {
					HashMap<String, String> map=split(eachline);
					String precursorMZString = map.get("Parent");
					if (precursorMZString!=null) {
						precursorMZ=Double.parseDouble(precursorMZString);
					}
					String scoreString=map.get("Unassigned");
					
					//Issue 90: If the score is missing, just write "0".
					//This is the same value used when importing a .blib generated in skyline.
					if (scoreString==null) {
						scoreString=map.get("Prob");
						if (scoreString != null) {
							score=Float.parseFloat(scoreString);
						} else {
							score = 0.0f; //Issue 90: force score of 0 if the score is missing.
						}
						
					} else {
						score=1.0f-Float.parseFloat(scoreString);
					}
					
					if (keepAccessions) {
						accession=map.get("Protein");
						if (accession!=null) {
							int stop=accession.indexOf(' ');
							if (stop>0) {
								accession=accession.substring(0, stop);
							}
						}
					}
					
					if (fullname==null) {
						fullname=map.get("Fullname");
						if (fullname == null) {
							fullname = altName; //Issue 90: If no other name information is included, fall back to sequence and charge.
						}
					}
					
					
					//Issue 90: sptext files do not always contain periods
					String sequence;
					String modLess = getModlessSequence(fullname);
					if (modLess.contains(".")) {
						sequence = modLess.substring(modLess.indexOf('.')+1, modLess.lastIndexOf('.'));
					} else {
						sequence = modLess;
					}
					
					sequence=PeptideUtils.getPeptideSeq(sequence);
					
					//Issue 90: The charge may be stored as a separate field in the comment,
					//or attached to the peptide sequence itself e.g. PEPTIDER+2
					if (map.containsKey("Charge")) {
						precursorCharge = Byte.parseByte(map.get("Charge"));
					} else if (precursorCharge!=(byte)0) {
						String substring = fullname.substring(fullname.lastIndexOf('/')+1);
						if (substring.indexOf(' ')>0) substring=substring.substring(0, substring.indexOf(' '));
						precursorCharge=Byte.parseByte(substring);
					}
					
					String mods=map.get("Mods");
					StringTokenizer st=new StringTokenizer(mods, "/()");
					int modCount=Integer.parseInt(st.nextToken());
					if (modCount>0) {
						TIntDoubleHashMap modMap=new TIntDoubleHashMap();
						while (st.hasMoreTokens()) {
							String modString=st.nextToken();
							StringTokenizer st2=new StringTokenizer(modString, ",");
							if (st2.countTokens()>2) {
								int index=Integer.parseInt(st2.nextToken());
								char aa=st2.nextToken().charAt(0);
								String mod=st2.nextToken().trim();
								double mass=getMass(aa, mod);
								
								if (modMap.contains(index)) {
									// shouldn't happen, but just in case
									modMap.put(index, modMap.get(index)+mass);
								} else {
									modMap.put(index, mass);
								}
							}
						}
						
						StringBuilder sb=new StringBuilder();
						for (int i=0; i<sequence.length(); i++) {
							sb.append(sequence.charAt(i));
							if (modMap.contains(i)) {
								sb.append('[');
								double mass=modMap.get(i);
								if (mass>=0) {
									sb.append('+');
								}
								sb.append(mass);
								sb.append(']');
							}
						}
						peptideModSeq=sb.toString();
						
					} else {
						peptideModSeq=sequence;
					}

					// uses peptideModSeq so needs to be last
					String rtString=map.get("RetentionTime");
					if (rtString!=null) {
						if (rtString.indexOf(',')>0) {
							rtString=rtString.substring(0, rtString.indexOf(','));
						}
						retentionTime=Float.parseFloat(rtString);
					} else if (retentionTime==0.0f) {
						retentionTime=(float)SSRCalc.getHydrophobicity(peptideModSeq);
					}
				}
			}
			if (peaks.size()>0) {
				Pair<double[], float[]> peakArrays=Peak.toArrays(peaks);
				HashSet<String> accessions=new HashSet<String>();
				if (accession!=null) accessions.add(accession);
				
				if (retentionTime==0.0f) Logger.errorLine("MSP parsing error: expected to find retention time for "+peptideModSeq+" but it was missing."); 
				if (precursorMZ==0.0) Logger.errorLine("MSP parsing error: expected to find precursor M/Z for "+peptideModSeq+" but it was missing."); 
				if (precursorCharge==(byte)0) Logger.errorLine("MSP parsing error: expected to find charge for "+peptideModSeq+" but it was missing."); 
				LibraryEntry entry=new LibraryEntry(fileName, accessions, precursorMZ, precursorCharge, peptideModSeq, 1, retentionTime, score, peakArrays.x, peakArrays.y, aaConstants);
				retentionTime=0.0f;
				precursorMZ=0.0;
				precursorCharge=0;
				
				entryList.add(entry);
			}
			return entryList;

		} finally {
			if (in!=null) {
				try {
					in.close();
				} catch (IOException ioe) {
					throw(ioe);
				}
			}
		}
	}

	/**
	 * Issue 90: Remove modifications given in brackets from a string,
	 * without examining contents of brackets.
	 * If no bracketed expressions exist, return original string
	 * <p>
	 * eg
	 * PEPT[80.0]IDER --> PEPTIDER
	 * [A]PEPT[80.0]IDER[B] --> PEPTIDER
	 * R.AAAAAAAAAAAAAAAGAGAGAK.Q --> R.AAAAAAAAAAAAAAAGAGAGAK.Q
	 * <p>
	 * @param sequence
	 * @return
	 */
	public static String getModlessSequence(String sequence) {
		
		StringBuilder sb = new StringBuilder();
		boolean isInMod = false;
		
		for (int i = 0; i < sequence.length(); i++) {
			
			char c = sequence.charAt(i);
			
			if (isInMod) {
				if (c == ']') {
					isInMod = false;
				}
				continue;
			} 
			
			if(!isInMod && c == '[') {
				isInMod = true;
				continue;
			} 
			
			sb.append(String.valueOf(c));
		}
		
		return sb.toString();
	}
	
	/**
	 *  inner quoted splitting because MSP format is insane
	 * @param s
	 * @return
	 */
	public static HashMap<String, String> split(String s) {
		char[] ca=s.toCharArray();
		ArrayList<String> strings=new ArrayList<String>();
		StringBuilder sb=new StringBuilder();
		boolean insideQuotes=false;
		for (int i=0; i<ca.length; i++) {
			if (ca[i]=='\"') {
				insideQuotes=!insideQuotes;
			} else if (ca[i]==' '&&!insideQuotes) {
				strings.add(sb.toString());
				sb.setLength(0);
			} else {
				sb.append(ca[i]);
			}
		}
		if (sb.length()>0) {
			strings.add(sb.toString());
		}
		
		HashMap<String, String> map=new HashMap<String, String>();
		for (String string : strings) {
			int equals=string.indexOf('=');
			if (equals>0) {
				String key=string.substring(0, equals);
				String value=string.substring(equals+1);
				map.put(key, value);
			} else {
				map.put(string, null);
			}
		}
		return map;
	}
	
	/**
	 * TODO should consider getting these out of a Unimod database. These are consistent with the NIST libraries as of January 2016
	 * @param aa
	 * @param mod
	 * @return
	 */
	public static double getMass(char aa, String mod) {
		if (aa=='C'&&"CAM".equals(mod)) {
			return 57.0214635;
		} else if (aa=='M'&&"Oxidation".equalsIgnoreCase(mod)) {
			return 	15.994915;
		} else if (aa=='C'&&"Methylthio".equalsIgnoreCase(mod)) {
			return 	45.987721;
		} else if (aa=='E'&&"Pyro_glu".equalsIgnoreCase(mod)) { // note "-" vs "_" (Unimod is crazy)
			return -18.010565;
		} else if (aa=='Q'&&"Pyro-glu".equalsIgnoreCase(mod)) { // note "-" vs "_" (Unimod is crazy)
			return -17.026549;
		} else if (aa=='E'&&"Glu->pyro-Glu".equalsIgnoreCase(mod)) {
			return -18.010565;
		} else if (aa=='Q'&&"Gln->pyro-Glu".equalsIgnoreCase(mod)) {
			return -17.026549;
		} else if (aa=='C'&&"Carbamidomethyl".equalsIgnoreCase(mod)) {
			return 	57.0214635;
		} else if (aa=='C'&&("Pyro-carbamidomethyl".equalsIgnoreCase(mod)||"Pyro-cmC".equalsIgnoreCase(mod))) {
			return 39.994915; // +57,-17
		} else if ("Acetyl".equalsIgnoreCase(mod)) {
			return 42.010565;
		} else if ("Deamidated".equalsIgnoreCase(mod)) {
			return 0.984016;
		}
		throw new EncyclopediaException("Unexpected modification ["+mod+"] on ["+aa+"]");
	}
}
