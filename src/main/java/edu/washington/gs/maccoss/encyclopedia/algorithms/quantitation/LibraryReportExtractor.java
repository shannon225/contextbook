package edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.TreeMap;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideReportData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroupInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroupQuantifier;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.QuantitativeDIAData;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;

public class LibraryReportExtractor {

	public static void main(String[] args) throws IOException, SQLException, DataFormatException {
		File file=new File("/Volumes/WorkingDisk/yeast_curve_wide.elib");

		LibraryFile library=new LibraryFile();
		library.openFile(file);
		
		LibraryReportExtractor.extractMatrix(library, true);
	}

	public static File extractMatrix(LibraryFile library, boolean normalizeByTIC) throws IOException, SQLException, DataFormatException {
		ArrayList<ProteinGroupInterface> proteins=library.getProteinGroups();
		return extractMatrix(library, proteins, normalizeByTIC, Optional.ofNullable(null));
	}
	public static File extractMatrix(LibraryFile library, ArrayList<ProteinGroupInterface> proteins, boolean normalizeByTIC) throws IOException, SQLException, DataFormatException {
		return extractMatrix(library, proteins, normalizeByTIC, Optional.ofNullable(null));
	}
	public static File extractMatrix(LibraryFile library, ArrayList<ProteinGroupInterface> proteins, boolean normalizeByTIC, Optional<CoefficientOfVariationCalculator> cvCalculator) throws IOException, SQLException, DataFormatException {
		File stubFile=library.getFile();
		if (stubFile==null) {
			throw new EncyclopediaException("Please save .ELIB before trying to read matrix data from it!");
		}
		String tag;
		if (normalizeByTIC) {
			tag="";
		} else {
			tag="_unormalized";
		}
		File diffactoFile=new File(stubFile.getParentFile(), stubFile.getName()+tag+".diffacto.csv");
		File diffactoSamplesFile=new File(stubFile.getParentFile(), stubFile.getName()+tag+".diffacto.samples.lst");
		File peptideReportFile=new File(stubFile.getParentFile(), stubFile.getName()+tag+".peptides.txt");
		File proteinReportFile=new File(stubFile.getParentFile(), stubFile.getName()+tag+".proteins.txt");
		
		Connection c=library.getConnection();
		try {
			Statement s=c.createStatement();
			PrintWriter diffactoWriter=null;
			PrintWriter diffactoSamplesWriter=null;
			PrintWriter peptideWriter=null;
			PrintWriter proteinWriter=null;
			try {
				ArrayList<String> sourceFiles=new ArrayList<String>();
				if (cvCalculator.isPresent()) {
					Logger.logLine("Using pre-selected source files...");
					sourceFiles=cvCalculator.get().getSortedSampleNames();
				} else {
					Logger.logLine("Getting source files...");
					ResultSet rs=s.executeQuery("select distinct SourceFile from peptidequants");
					while (rs.next()) {
						sourceFiles.add(rs.getString(1));
					}
					rs.close();

					Collections.sort(sourceFiles);
				}
				TObjectIntHashMap<String> indexByFile=new TObjectIntHashMap<>();
				for (int i=0; i<sourceFiles.size(); i++) {
					indexByFile.put(sourceFiles.get(i), i);
				}
				
				ArrayList<ProteinGroupQuantifier> proteinQuantifiers=new ArrayList<ProteinGroupQuantifier>();
				for (int i=0; i<sourceFiles.size(); i++) {
					proteinQuantifiers.add(new ProteinGroupQuantifier(proteins));
				}
				
				diffactoSamplesWriter=new PrintWriter(diffactoSamplesFile, "UTF-8");
				
				diffactoWriter=new PrintWriter(diffactoFile, "UTF-8");
				diffactoWriter.print("sequences,Protein");
				
				peptideWriter=new PrintWriter(peptideReportFile, "UTF-8");
				peptideWriter.print("Peptide\tProtein\tnumFragments");
				
				proteinWriter=new PrintWriter(proteinReportFile, "UTF-8");
				proteinWriter.print("Protein\tNumPeptides\tPeptideSequences");
				
				float averageTIC=0.0f;
				TObjectFloatHashMap<String> ticBySourceFileMap=new TObjectFloatHashMap<String>();
				int sampleGroup=0;
				for (String sourceFile : sourceFiles) {
					sampleGroup++;
					if (normalizeByTIC) {
						float tic=library.getTIC(sourceFile);
						ticBySourceFileMap.put(sourceFile, tic);
						averageTIC+=tic;
					}
					diffactoWriter.print(","+sourceFile);
					peptideWriter.print("\t"+sourceFile);
					proteinWriter.print("\t"+sourceFile);
					diffactoSamplesWriter.print(sourceFile+"\tS"+sampleGroup);
					diffactoSamplesWriter.println();
				}
				if (sourceFiles.size()>0) {
					averageTIC=averageTIC/sourceFiles.size();
				}
				
				diffactoWriter.println();
				peptideWriter.println();
				proteinWriter.println();
				Logger.logLine("Found "+sourceFiles.size()+" data files");
				
				HashMap<String, int[]> numFragmentsByPeptideModSeq=new HashMap<String, int[]>();
				TreeMap<String, Pair<String, float[]>> intensitiesByPeptideModSeq=new TreeMap<String, Pair<String, float[]>>();
				ResultSet rs = s.executeQuery("select " +
						"pep.PrecursorCharge, " +
						"pep.PeptideModSeq, " +
						"pep.SourceFile, " +
						"pep.TotalIntensity, " +
						"pep.NumberOfQuantIons, " +
						"group_concat(p.ProteinAccession, '" + PSMData.ACCESSION_TOKEN + "') as ProteinAccessions " +
						"from " +
						"peptidequants pep " +
						"join peptidescores s using (peptidemodseq, precursorcharge) " + // outer join to scores table means we'll skip quant rows from unscored charge states
						"left join peptidetoprotein p " +
						"where " +
						"pep.PeptideSeq = p.PeptideSeq " +
						"group by pep.rowid;"
				);
				int count=0;
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
					
					int index=indexByFile.get(sourceFile);
					if (index<0) {
						Logger.errorLine("Can't find ["+sourceFile+"]!");
						Logger.errorLine("Keys: {");
						for (String name : sourceFiles) {
							Logger.errorLine("    "+name);
						}
						Logger.errorLine("}");
						throw new EncyclopediaException("Unexpected sample: "+sourceFile);
					}

					// FIXME NEED TO NORMALIZE BY TIC
					float normalizedIntensity;
					if (cvCalculator.isPresent()) {
						normalizedIntensity=totalIntensity*cvCalculator.get().getReplicateNormalizationFactor(sourceFile, ticBySourceFileMap);
					} else {
						float tic=ticBySourceFileMap.get(sourceFile);
						if (tic>0.0f) {
							normalizedIntensity=totalIntensity/tic*averageTIC;
						} else {
							normalizedIntensity=totalIntensity;
						}
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
				
				int totalAdded=0;
				HashMap<String, ArrayList<String>> peptidesInProtein=new HashMap<>();
				HashSet<String> badCVPeptides=new HashSet<>();
				HashSet<String> badCompletenessPeptides=new HashSet<>();
				TObjectFloatHashMap<String> cvs=new TObjectFloatHashMap<>();
				for (Entry<String, Pair<String, float[]>> entry : intensitiesByPeptideModSeq.entrySet()) {
					String peptideModSeq=entry.getKey();
					String proteinToken=entry.getValue().x;
					HashSet<String> accessions=PSMData.stringToAccessions(proteinToken);
					float[] intensitiesArray=entry.getValue().y;
					
					if (cvCalculator.isPresent()) {
						Pair<Float, Boolean> pair=cvCalculator.get().getCV(sourceFiles, intensitiesArray);
						float cv=pair.x;
						boolean atLeastSampleFullyMeasured=pair.y;
						
						cvs.put(peptideModSeq, cv);
						if (Float.isNaN(cv)||cv>cvCalculator.get().getMaximumAcceptedCV()) {
							badCVPeptides.add(peptideModSeq);
							continue;
						};
						
						if (!atLeastSampleFullyMeasured) {
							badCompletenessPeptides.add(peptideModSeq);
							continue;
						}
					}
					
					boolean anyAdded=false;
					for (int index=0; index<intensitiesArray.length; index++) {
						boolean added=proteinQuantifiers.get(index).addIntensity(accessions, intensitiesArray[index]);
						if (added) {
							anyAdded=true;
							totalAdded++;
						}
					}
					
					if (anyAdded) {
						for (String accession : accessions) {
							ArrayList<String> peptides=peptidesInProtein.get(accession);
							if (peptides==null) {
								peptides=new ArrayList<>();
								peptidesInProtein.put(accession, peptides);
							}
							peptides.add(peptideModSeq);
						}
						
					}
				}
				
				if (cvCalculator.isPresent()) {
					Logger.logLine("Finished processing "+totalAdded+"/"+count+" measurements. Found "+intensitiesByPeptideModSeq.size()+" quantitative peptides (where "+badCVPeptides.size()+" were outside a CV of "+cvCalculator.get().getMaximumAcceptedCV()+" and "+badCompletenessPeptides.size()+" were incomplete). Writing reports...");

					File cvReportFile=new File(stubFile.getParentFile(), stubFile.getName()+".cvs.txt");
					PrintWriter cvWriter=null; 
					try {
						cvWriter=new PrintWriter(cvReportFile, "UTF-8");
						cvWriter.println("peptide\tcv");
						
						final PrintWriter finalCVWriter=cvWriter;
						cvs.forEachEntry(new TObjectFloatProcedure<String>() {
							@Override
							public boolean execute(String a, float b) {
								finalCVWriter.println(a+"\t"+b);
								return true;
							}
						});
					} finally {
						if (cvWriter!=null) {
							cvWriter.close();
						}
					}
				} else {
					Logger.logLine("Finished processing "+count+" records, found "+intensitiesByPeptideModSeq.size()+" quantitative unique peptides. Writing reports...");
				}
				
				int numberInconsistentFragments=0;
				for (Entry<String, Pair<String, float[]>> entry : intensitiesByPeptideModSeq.entrySet()) {
					String peptideModSeq=entry.getKey();
					if (badCompletenessPeptides.size()>0&&badCompletenessPeptides.contains(peptideModSeq)) {
						continue;
					} else if (badCVPeptides.size()>0&&badCVPeptides.contains(peptideModSeq)) {
						continue;
					}
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
					
					for (float f : array) {
						peptideWriter.print("\t");
						peptideWriter.print(f);
					}
					peptideWriter.println();
					
					// DIFFACTO
					diffactoWriter.print("//");
					diffactoWriter.print(peptideModSeq);
					diffactoWriter.print(",");
					ArrayList<String> accessions=new ArrayList<>(PSMData.stringToAccessions(pair.x));
					Collections.sort(accessions);
					diffactoWriter.print(General.toString(accessions.toArray(), ";"));

					for (float f : array) {
						diffactoWriter.print(",");
						diffactoWriter.print(f);
					}
					diffactoWriter.println();
				}
				if (numberInconsistentFragments>0) {
					Logger.errorLine("Inconsistent number of fragments in "+numberInconsistentFragments+" of "+intensitiesByPeptideModSeq.size()+" peptides");
				}
				Logger.logLine("Finished writing peptide report for "+totalAdded+" unique peptides!");
				
				int numberProteinsKept=0;
				for (ProteinGroupInterface protein : proteins) {
					if (protein.isDecoy()) continue;
					
					HashSet<String> peptides=new HashSet<>();
					for (String accession : protein.getEquivalentAccessions()) {
						ArrayList<String> thisAccessionsPeptides=peptidesInProtein.get(accession);
						if (thisAccessionsPeptides!=null) {
							peptides.addAll(thisAccessionsPeptides);
						}
					}
					ArrayList<String> sortedPeptides=new ArrayList<>(peptides);
					Collections.sort(sortedPeptides);
					
					StringBuilder sb=new StringBuilder(protein.toString());
					sb.append("\t"+sortedPeptides.size());
					boolean first=true;
					for (String peptide : sortedPeptides) {
						if (first) {
							sb.append("\t");
							first=false;
						} else {
							sb.append(";");
						}
						sb.append(peptide);
					}
					float totalIntensity=0.0f;
					for (ProteinGroupQuantifier proteinQuantifier : proteinQuantifiers) {
						float intensity=proteinQuantifier.getIntensity(protein);
						totalIntensity+=intensity;
						sb.append("\t"+intensity);
					}
					
					if (totalIntensity>0.0f) {
						numberProteinsKept++;
						proteinWriter.println(sb.toString());
					}
				}
				Logger.logLine("Finished writing protein report for "+numberProteinsKept+" protein groups!");
				
				rs.close();
			} finally {
				s.close();
				if (diffactoSamplesWriter!=null) diffactoSamplesWriter.close();
				if (diffactoWriter!=null) diffactoWriter.close();
				if (peptideWriter!=null) peptideWriter.close();
				if (proteinWriter!=null) proteinWriter.close();
			}
		} finally {
			c.close();
		}
		return proteinReportFile;
	}

	public static Pair<ArrayList<String>, ArrayList<PeptideReportData>> extractMatrix(LibraryFile library, AminoAcidConstants aaConstants) throws IOException, SQLException, DataFormatException {
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
				
				// FIXME need to add TIC normalization in
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
				
				rs = s.executeQuery("select " +
						"pep.PrecursorCharge, " +
						"pep.PeptideModSeq, " +
						"pep.SourceFile, " +
						"pep.RTInSecondsCenter, pep.RTInSecondsStart, pep.RTInSecondsStop, " +
						"group_concat(p.ProteinAccession, '" + PSMData.ACCESSION_TOKEN + "') as ProteinAccessions, " +
						"pep.QuantIonMassLength, pep.QuantIonMassArray, "+
						"pep.QuantIonIntensityLength, pep.QuantIonIntensityArray "+
						"from " +
						"peptidequants pep " +
						"join peptidescores s using (peptidemodseq, precursorcharge) " + // outer join to scores table means we'll skip quant rows from unscored charge states
						"left join peptidetoprotein p " +
						"where " +
						"pep.PeptideSeq = p.PeptideSeq " +
						"group by pep.rowid;"
				);
				
				int count=0;
				while (rs.next()) {
					count++;
					if (count%10000==0) {
						Logger.logLine(count+" records processed...");
					}
					byte precursorCharge=rs.getByte(1);
					String peptideModSeq=rs.getString(2);
					String sourceFile=rs.getString(3);
					float scanStartTime=rs.getFloat(4);
					float rtStart=rs.getFloat(5);
					float rtStop=rs.getFloat(6);
					String accessions=rs.getString(7);
					
					int index=Collections.binarySearch(sourceFiles, sourceFile);
					if (index<0) throw new EncyclopediaException("Unexpected sample: "+sourceFile);

					PeptideReportData data=intensitiesByPeptideModSeq.get(peptideModSeq);
					if (data==null) {
						data=new PeptideReportData(peptideModSeq, precursorCharge, accessions, aaConstants);
						intensitiesByPeptideModSeq.put(peptideModSeq, data);
					}

					Range rtScanRange=new Range(rtStart, rtStop);
					int quantIonMassesLength=rs.getInt(8);
					double[] quantIonMasses;
					if (quantIonMassesLength>0) {
						quantIonMasses=ByteConverter.toDoubleArray(CompressionUtils.decompress(rs.getBytes(9), quantIonMassesLength));
					} else {
						quantIonMasses=new double[] {};
					}
					int quantIonIntensityLength=rs.getInt(10);
					float[] quantIonIntensities;
					if (quantIonIntensityLength>0) {
						quantIonIntensities=ByteConverter.toFloatArray(CompressionUtils.decompress(rs.getBytes(11), quantIonIntensityLength));
					} else {
						quantIonIntensities=new float[] {};
					}
					QuantitativeDIAData quantData=new QuantitativeDIAData(peptideModSeq, precursorCharge, scanStartTime, rtScanRange, quantIonMasses, quantIonIntensities, aaConstants);
					data.addQuantitativeDIAData(sourceFile, quantData);
				}
				Logger.logLine("Finished processing "+count+" records, found "+intensitiesByPeptideModSeq.size()+" quantitative unique peptides. Writing reports...");
				
				ArrayList<PeptideReportData> reportData=new ArrayList<PeptideReportData>(intensitiesByPeptideModSeq.values());
				
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
