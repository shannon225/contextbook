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
import edu.washington.gs.maccoss.encyclopedia.datastructures.*;
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
		
		convertMSP(mspFile, fastaFile, libraryFile);
	}
	public static void convertMSP(File mspFile, File fastaFile) throws IOException, SQLException {
		String absolutePath=mspFile.getAbsolutePath();
		File libraryFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+LibraryFile.DLIB);
		convertMSP(mspFile, fastaFile, libraryFile);
	}

	public static void convertMSP(File mspFile, File fastaFile, File libraryFile) throws IOException, SQLException {
		Logger.logLine("Reading MSP file "+mspFile.getName());
		ArrayList<LibraryEntry> entries=readMSP(mspFile, false);

		Logger.logLine("Reading Fasta file "+fastaFile.getName());
		ArrayList<FastaEntryInterface> proteins=FastaReader.readFasta(fastaFile);

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
		int batchSize=entries.size()/10;
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
	
	public static ArrayList<LibraryEntry> readMSP(File f, boolean keepAccessions) {
		BufferedReader in=null;
		ArrayList<LibraryEntry> entryList=new ArrayList<LibraryEntry>();
		try {
			in=new BufferedReader(new FileReader(f));
			return readMSP(in, f.getName(), keepAccessions);

		} catch (IOException ioe) {
			Logger.errorLine("I/O Error found reading NIST MSP Library ["+f.getAbsolutePath()+"]");
			Logger.errorException(ioe);
			return entryList;
		} finally {
			if (in!=null) {
				try {
					in.close();
				} catch (IOException ioe) {
					ioe.printStackTrace();
				}
			}
		}
	}
	
	public static ArrayList<LibraryEntry> readMSP(String s, String fileName, boolean keepAccessions) {
		return readMSP(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(s.getBytes(StandardCharsets.UTF_8)))), fileName, keepAccessions);
	}
	
	public static ArrayList<LibraryEntry> readMSP(InputStream s, String fileName, boolean keepAccessions) {
		return readMSP(new BufferedReader(new InputStreamReader(s)), fileName, keepAccessions);
	}
	
	public static ArrayList<LibraryEntry> readMSP(BufferedReader in, String fileName, boolean keepAccessions) {
		//TODO: take in parameters
		final AminoAcidConstants aaConstants = new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());

		ArrayList<LibraryEntry> entryList=new ArrayList<LibraryEntry>();
		String eachline=null;
		try {
			
			String peptideModSeq=null;
			String fullname=null;
			String accession=null;
			byte precursorCharge=0;
			double precursorMZ=0.0;
			float retentionTime=0.0f;
			float score=0.0f;
			ArrayList<Peak> peaks=new ArrayList<Peak>();

			OUTERLOOP: while ((eachline=in.readLine())!=null) {
				if (eachline.trim().length()==0) {
					continue OUTERLOOP;
				}
				if (eachline.startsWith("Name: ")) {
					if (peaks.size()>0) {
						Pair<double[], float[]> peakArrays=Peak.toArrays(peaks);
						HashSet<String> accessions=new HashSet<String>();
						if (accession!=null) accessions.add(accession);
						LibraryEntry entry=new LibraryEntry(fileName, accessions, precursorMZ, precursorCharge, peptideModSeq, 1, retentionTime, score, peakArrays.x, peakArrays.y, aaConstants);
						entryList.add(entry);
						peaks.clear();
						fullname=null;
						
					}
				} else if (eachline.startsWith("Num peaks: ")||eachline.startsWith("NumPeaks: ")) {
					INNERLOOP: while ((eachline=in.readLine())!=null) {
						if (eachline.trim().length()==0) {
							break INNERLOOP;
						}
						StringTokenizer st=new StringTokenizer(eachline);
						double mass=Double.parseDouble(st.nextToken());
						float intensity=Float.parseFloat(st.nextToken());
						peaks.add(new Peak(mass, intensity));
					}
				} else if (eachline.startsWith("FullName: ")) {
					fullname=eachline.substring(10);
				} else if (eachline.startsWith("Comment: ")) {
					HashMap<String, String> map=split(eachline);
					precursorMZ=Double.parseDouble(map.get("Parent"));
					String scoreString=map.get("Unassigned");
					
					if (scoreString==null) {
						scoreString=map.get("Prob");
						score=Float.parseFloat(scoreString);
						
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

					String rtString=map.get("RetentionTime");
					if (rtString!=null) {
						if (rtString.indexOf(',')>0) {
							rtString=rtString.substring(0, rtString.indexOf(','));
						}
						retentionTime=Float.parseFloat(rtString);
					} else {
						retentionTime=(float)SSRCalc.getHydrophobicity(peptideModSeq);
					}
					
					if (fullname==null) {
						fullname=map.get("Fullname");
					}
					String sequence=fullname.substring(fullname.indexOf('.')+1, fullname.lastIndexOf('.'));
					sequence=PeptideUtils.getPeptideSeq(sequence);
					precursorCharge=Byte.parseByte(fullname.substring(fullname.lastIndexOf('/')+1));
					
					String mods=map.get("Mods");
					StringTokenizer st=new StringTokenizer(mods, "/");
					int modCount=Integer.parseInt(st.nextToken());
					if (modCount>0) {
						TIntDoubleHashMap modMap=new TIntDoubleHashMap();
						while (st.hasMoreTokens()) {
							StringTokenizer st2=new StringTokenizer(st.nextToken(), ",");
							int index=Integer.parseInt(st2.nextToken());
							char aa=st2.nextToken().charAt(0);
							String mod=st2.nextToken();
							double mass=getMass(aa, mod);
							
							if (modMap.contains(index)) {
								// shouldn't happen, but just in case
								modMap.put(index, modMap.get(index)+mass);
							} else {
								modMap.put(index, mass);
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
				}
			}
			if (peaks.size()>0) {
				Pair<double[], float[]> peakArrays=Peak.toArrays(peaks);
				HashSet<String> accessions=new HashSet<String>();
				if (accession!=null) accessions.add(accession);
				LibraryEntry entry=new LibraryEntry(fileName, accessions, precursorMZ, precursorCharge, peptideModSeq, 1, retentionTime, score, peakArrays.x, peakArrays.y, aaConstants);
				entryList.add(entry);
			}
			return entryList;

		} catch (IOException ioe) {
			Logger.errorLine("I/O Error found reading NIST MSP Library ["+fileName+"]");
			Logger.errorException(ioe);
			return entryList;
		} catch (Exception e) {
			Logger.errorLine("I/O Error found reading NIST MSP Library ["+fileName+"], parsing ["+eachline+"]");
			Logger.errorException(e);
			return entryList;
		} finally {
			if (in!=null) {
				try {
					in.close();
				} catch (IOException ioe) {
					ioe.printStackTrace();
				}
			}
		}
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
		} else if (aa=='C'&&"Pyro-carbamidomethyl".equalsIgnoreCase(mod)) {
			return 39.994915; // +57,-17
		} else if ("Acetyl".equalsIgnoreCase(mod)) {
			return 42.010565;
		}
		throw new EncyclopediaException("Unexpected modification ["+mod+"] on ["+aa+"]");
	}
}
