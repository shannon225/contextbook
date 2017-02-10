package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroup;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroupQuantifier;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class LibraryReportExtractor {
	public static void extractMatrix(LibraryFile library, ArrayList<ProteinGroup> proteins) throws IOException, SQLException, DataFormatException {
		File stubFile=library.getFile();
		if (stubFile==null) {
			throw new EncyclopediaException("Please save .ELIB before trying to read matrix data from it!");
		}
		File peptideReportFile=new File(stubFile.getParentFile(), stubFile.getName()+".peptides.txt");
		File proteinReportFile=new File(stubFile.getParentFile(), stubFile.getName()+".proteins.txt");
		
		Connection c=library.getConnection();
		try {
			Statement s=c.createStatement();
			PrintWriter peptideWriter=null;
			PrintWriter proteinWriter=null;
			try {
				ArrayList<String> sourceFiles=new ArrayList<String>();
				
				Logger.logLine("Getting source files...");
				ResultSet rs=s.executeQuery("select distinct SourceFile from peptidequants");
				while (rs.next()) {
					sourceFiles.add(rs.getString(1));
				}
				rs.close();
				
				Collections.sort(sourceFiles);
				
				ArrayList<ProteinGroupQuantifier> proteinQuantifiers=new ArrayList<ProteinGroupQuantifier>();
				for (int i=0; i<sourceFiles.size(); i++) {
					proteinQuantifiers.add(new ProteinGroupQuantifier(proteins));
				}
				
				peptideWriter=new PrintWriter(peptideReportFile, "UTF-8");
				peptideWriter.print("Peptide\tProtein\tnumFragments");
				
				proteinWriter=new PrintWriter(proteinReportFile, "UTF-8");
				proteinWriter.print("Protein\tnumEquivalentAccessions");
				
				float averageTIC=0.0f;
				TObjectFloatHashMap<String> ticBySourceFileMap=new TObjectFloatHashMap<String>();
				for (String sourceFile : sourceFiles) {
					float tic=library.getTIC(sourceFile);
					ticBySourceFileMap.put(sourceFile, tic);
					averageTIC+=tic;
					
					peptideWriter.print("\t"+sourceFile);
					proteinWriter.print("\t"+sourceFile);
					proteinWriter.print("\tnum_"+sourceFile);
				}
				averageTIC=averageTIC/sourceFiles.size();
				
				peptideWriter.println();
				proteinWriter.println();
				Logger.logLine("Found "+sourceFiles.size()+" data files");
				
				HashMap<String, int[]> numFragmentsByPeptideModSeq=new HashMap<String, int[]>();
				TreeMap<String, Pair<String, float[]>> intensitiesByPeptideModSeq=new TreeMap<String, Pair<String, float[]>>();
				rs=s.executeQuery("select pep.PrecursorCharge, pep.PeptideModSeq, pep.SourceFile, pep.TotalIntensity, pep.NumberOfQuantIons, pro.ProteinAccessions from peptidequants pep, proteins pro where pep.PeptideSeq = pro.PeptideSeq");
				int count=0;
				int totalAdded=0;
				while (rs.next()) {
					count++;
					if (count%10000==0) {
						Logger.logLine(count+" records processed...");
					}
					//byte precursorCharge=(byte)rs.getInt(1);
					String peptideModSeq=rs.getString(2);
					String sourceFile=rs.getString(3);
					float totalIntensity=rs.getFloat(4);
					int numberOfQuantIons=rs.getInt(5);
					String proteinToken=rs.getString(6);
					HashSet<String> accessions=PSMData.stringToAccessions(proteinToken);
					
					int index=Collections.binarySearch(sourceFiles, sourceFile);
					if (index<0) throw new EncyclopediaException("Unexpected sample: "+sourceFile);

					float tic=ticBySourceFileMap.get(sourceFile);
					float normalizedIntensity;
					if (tic>0.0f) {
						normalizedIntensity=totalIntensity/tic*averageTIC;
					} else {
						normalizedIntensity=totalIntensity;
					}
					
					boolean added=proteinQuantifiers.get(index).addIntensity(accessions, normalizedIntensity);
					if (added) {
						totalAdded++;
					}
					
					Pair<String, float[]> pair=intensitiesByPeptideModSeq.get(peptideModSeq);
					int[] numFragmentsArray=numFragmentsByPeptideModSeq.get(peptideModSeq);
					float[] intensitiesArray;
					if (pair==null) {
						intensitiesArray=new float[sourceFiles.size()];
						intensitiesByPeptideModSeq.put(peptideModSeq, new Pair<String, float[]>(proteinToken, intensitiesArray));
						numFragmentsArray=new int[sourceFiles.size()];
						numFragmentsByPeptideModSeq.put(peptideModSeq, numFragmentsArray);
					} else {
						intensitiesArray=pair.y;
					}
					intensitiesArray[index]+=normalizedIntensity; // sums charge states together
					numFragmentsArray[index]=numberOfQuantIons;
				}
				Logger.logLine("Finished processing "+count+" records, found "+totalAdded+" quantitative unique peptides. Writing reports...");
				
				int numberInconsistentFragments=0;
				for (Entry<String, Pair<String, float[]>> entry : intensitiesByPeptideModSeq.entrySet()) {
					String peptideModSeq=entry.getKey();
					Pair<String, float[]> pair=entry.getValue();
					peptideWriter.print(peptideModSeq);
					peptideWriter.print("\t");
					peptideWriter.print(pair.x);
					int[] numFragments=numFragmentsByPeptideModSeq.get(peptideModSeq);
					
					int maxNumFragments=General.max(numFragments);
					int minNumFragments=General.min(numFragments);
					if (minNumFragments!=maxNumFragments) {
						numberInconsistentFragments++;
					}
					peptideWriter.print("\t");
					peptideWriter.print(maxNumFragments);
					
					float[] array=pair.y;
					
					/*if (array.length>1) {
						float mean=General.mean(array);
						float stdev=General.stdev(array);
						for (int i=0; i<array.length; i++) {
							if (stdev==0) {
								array[i]=0.0f;
							} else {
								array[i]=(array[i]-mean)/stdev;
							}
						}
					}*/
					
					for (float f : array) {
						peptideWriter.print("\t");
						peptideWriter.print(f);
					}
					peptideWriter.println();
				}
				if (numberInconsistentFragments>0) {
					Logger.errorLine("Inconsistent number of fragments in "+numberInconsistentFragments+" of "+intensitiesByPeptideModSeq.size()+" peptides");
				}
				Logger.logLine("Finished writing peptide report!");
				
				for (ProteinGroup protein : proteins) {
					//System.out.println(protein.getEquivalentAccessions().size()+"\t"+protein.getNspScore()+"\t"+protein.toString());
					proteinWriter.print(protein.toString());
					proteinWriter.print("\t"+protein.getEquivalentAccessions().size()); // numPeptides
					for (ProteinGroupQuantifier proteinQuantifier : proteinQuantifiers) {
						float intensity=proteinQuantifier.getIntensity(protein);
						proteinWriter.print("\t"+intensity);
						proteinWriter.print("\t"+proteinQuantifier.getNumberOfQuantitativePeptides(protein));
					}
					proteinWriter.println();
				}
				Logger.logLine("Finished writing protein report!");
				
				rs.close();
			} finally {
				s.close();
				if (peptideWriter!=null) peptideWriter.close();
				if (proteinWriter!=null) proteinWriter.close();
			}
		} finally {
			c.close();
		}
	}
	
	public static class PeptideReportData {
		private final String peptideModSeq;
		private final byte precursorCharge;
		private final String accessions;
		private final Range[] rtRanges;
		private final float[] totalIntensities;
		private final double[] targetFragmentMzs;
		private float avgRT=-1;
		public PeptideReportData(String peptideModSeq, byte precursorCharge, String accessions, double[] targetFragmentMzs, Range[] rtRanges, float[] totalIntensities) {
			this.peptideModSeq=peptideModSeq;
			this.precursorCharge=precursorCharge;
			this.accessions=accessions;
			this.rtRanges=rtRanges;
			this.totalIntensities=totalIntensities;
			this.targetFragmentMzs=targetFragmentMzs;
		}
		public String getPeptideModSeq() {
			return peptideModSeq;
		}
		public byte getPrecursorCharge() {
			return precursorCharge;
		}
		public String getAccessions() {
			return accessions;
		}
		public double[] getTargetFragmentMzs() {
			return targetFragmentMzs;
		}
		public Range[] getRTRanges() {
			return rtRanges;
		}
		public float[] getTotalIntensities() {
			return totalIntensities;
		}
		public float getAverageRetentionTime() {
			if (avgRT<0) {
				TFloatArrayList rtCenters=new TFloatArrayList();
				for (int i=0; i<rtRanges.length; i++) {
					rtCenters.add(rtRanges[i].getMiddle());
				}
				avgRT=General.mean(rtCenters.toArray());
			}
			return avgRT;
		}
	}

	public static Pair<ArrayList<String>, ArrayList<PeptideReportData>> extractMatrix(LibraryFile library) throws IOException, SQLException, DataFormatException {
		File stubFile=library.getFile();
		if (stubFile==null) {
			throw new EncyclopediaException("Please save .ELIB before trying to read matrix data from it!");
		}
		
		Connection c=library.getConnection();
		try {
			Statement s=c.createStatement();
			try {
				ArrayList<String> sourceFiles=new ArrayList<String>();
				
				Logger.logLine("Getting source files...");
				ResultSet rs=s.executeQuery("select distinct SourceFile from peptidequants");
				while (rs.next()) {
					sourceFiles.add(rs.getString(1));
				}
				rs.close();
				
				Collections.sort(sourceFiles);
				
				float averageTIC=0.0f;
				TObjectFloatHashMap<String> ticBySourceFileMap=new TObjectFloatHashMap<String>();
				for (String sourceFile : sourceFiles) {
					float tic=library.getTIC(sourceFile);
					ticBySourceFileMap.put(sourceFile, tic);
					averageTIC+=tic;
				}
				averageTIC=averageTIC/sourceFiles.size();
				
				Logger.logLine("Found "+sourceFiles.size()+" data files");
				
				TreeMap<String, PeptideReportData> intensitiesByPeptideModSeq=new TreeMap<String, PeptideReportData>();
				
				rs=s.executeQuery("select pep.PrecursorCharge, pep.PeptideModSeq, pep.SourceFile, pep.RTInSecondsStart, pep.RTInSecondsStop, pep.TotalIntensity, pro.ProteinAccessions, pep.QuantIonMassLength, pep.QuantIonMassArray from peptidequants pep, proteins pro where pep.PeptideSeq = pro.PeptideSeq");
				
				int count=0;
				int totalAdded=0;
				while (rs.next()) {
					count++;
					if (count%10000==0) {
						Logger.logLine(count+" records processed...");
					}
					byte precursorCharge=rs.getByte(1);
					String peptideModSeq=rs.getString(2);
					String sourceFile=rs.getString(3);
					float rtStart=rs.getFloat(4);
					float rtStop=rs.getFloat(5);
					float totalIntensity=rs.getFloat(6);
					String accessions=rs.getString(7);
					
					int index=Collections.binarySearch(sourceFiles, sourceFile);
					if (index<0) throw new EncyclopediaException("Unexpected sample: "+sourceFile);

					PeptideReportData data=intensitiesByPeptideModSeq.get(peptideModSeq);
					Range[] rtRange;
					float[] totalIntensities;
					if (data==null) {
						int quantIonMassesLength=rs.getInt(8);
						double[] quantIonMasses=ByteConverter.toDoubleArray(CompressionUtils.decompress(rs.getBytes(9), quantIonMassesLength));
						Arrays.sort(quantIonMasses);
						rtRange=new Range[sourceFiles.size()];
						totalIntensities=new float[sourceFiles.size()];
						intensitiesByPeptideModSeq.put(peptideModSeq, new PeptideReportData(peptideModSeq, precursorCharge, accessions, quantIonMasses, rtRange, totalIntensities));
					} else {
						rtRange=data.getRTRanges();
						totalIntensities=data.getTotalIntensities();
					}
					rtRange[index]=new Range(rtStart, rtStop);
					totalIntensities[index]=totalIntensity;
				}
				Logger.logLine("Finished processing "+count+" records, found "+totalAdded+" quantitative unique peptides. Writing reports...");
				
				ArrayList<PeptideReportData> reportData=new ArrayList<LibraryReportExtractor.PeptideReportData>(intensitiesByPeptideModSeq.values());
				
				Logger.logLine("Finished extracting peptide report!");
				
				rs.close();

				return new Pair<ArrayList<String>, ArrayList<PeptideReportData>>(sourceFiles, reportData);
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}
}
