package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

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
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroup;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroupQuantifier;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
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
				peptideWriter.print("Peptide");
				
				proteinWriter=new PrintWriter(proteinReportFile, "UTF-8");
				proteinWriter.println("Protein");
				
				float averageTIC=0.0f;
				TObjectFloatHashMap<String> ticBySourceFileMap=new TObjectFloatHashMap<String>();
				for (String sourceFile : sourceFiles) {
					float tic=library.getTIC(sourceFile);
					ticBySourceFileMap.put(sourceFile, tic);
					averageTIC+=tic;
					
					peptideWriter.print("\t"+sourceFile);
					proteinWriter.print("\t"+sourceFile);
				}
				averageTIC=averageTIC/sourceFiles.size();
				
				peptideWriter.println();
				proteinWriter.println();
				Logger.logLine("Found "+sourceFiles.size()+" data files");
				
				HashMap<String, float[]> intensitiesByPeptideModSeq=new HashMap<String, float[]>();
				rs=s.executeQuery("select pep.PrecursorCharge, pep.PeptideModSeq, pep.SourceFile, pep.TotalIntensity, pro.ProteinAccessions from peptidequants pep, proteins pro where pep.PeptideSeq = pro.PeptideSeq");
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
					String proteinToken=rs.getString(5);
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
					
					proteinQuantifiers.get(index).addIntensity(accessions, normalizedIntensity);
					
					float[] array=intensitiesByPeptideModSeq.get(peptideModSeq);
					if (array==null) {
						array=new float[sourceFiles.size()];
						intensitiesByPeptideModSeq.put(peptideModSeq, array);
					}
					array[index]+=normalizedIntensity; // sums charge states together
				}
				Logger.logLine("Finished processing "+count+" records, writing reports...");
				
				for (Entry<String, float[]> entry : intensitiesByPeptideModSeq.entrySet()) {
					peptideWriter.print(entry.getKey());
					float[] array=entry.getValue();
					
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
						peptideWriter.print("\t"+f);
					}
					peptideWriter.println();
				}
				Logger.logLine("Finished writing peptide report!");
				
				for (ProteinGroup protein : proteins) {
					peptideWriter.print(protein.toString());
					for (ProteinGroupQuantifier proteinQuantifier : proteinQuantifiers) {
						float intensity=proteinQuantifier.getIntensity(protein);
						proteinWriter.print("\t"+intensity);
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
}
