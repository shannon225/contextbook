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
import java.util.Map.Entry;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class LibraryReportExtractor {
	public static void main(String[] args) throws Exception {
		File libraryFile=new File("/Volumes/WorkingDisk/yeast_curve_wide.elib");
		File reportFile=new File("/Volumes/WorkingDisk/yeast_curve_wide.report.txt");
		
		LibraryInterface library=BlibToLibraryConverter.getFile(libraryFile);
		extractMatrix((LibraryFile)library, reportFile);
	}
	public static void extractMatrix(LibraryFile library, File reportFile) throws IOException, SQLException, DataFormatException {
		Connection c=library.getConnection();
		try {
			Statement s=c.createStatement();
			PrintWriter writer=null;
			try {
				ArrayList<String> sourceFiles=new ArrayList<String>();
				
				Logger.logLine("Getting source files...");
				ResultSet rs=s.executeQuery("select distinct SourceFile from peptidequants");
				while (rs.next()) {
					sourceFiles.add(rs.getString(1));
				}
				rs.close();
				
				Collections.sort(sourceFiles);
				
				writer=new PrintWriter(reportFile, "UTF-8");
				writer.print("Peptide");
				
				TObjectFloatHashMap<String> ticBySourceFileMap=new TObjectFloatHashMap<String>();
				for (String sourceFile : sourceFiles) {
					float tic=library.getTIC(sourceFile);
					ticBySourceFileMap.put(sourceFile, tic);
					
					writer.print("\t"+sourceFile);
				}
				writer.println();
				Logger.logLine("Found "+sourceFiles.size()+" data files");
				
				HashMap<String, float[]> intensitiesByPeptideModSeq=new HashMap<String, float[]>();
				rs=s.executeQuery("select PrecursorCharge, PeptideModSeq, SourceFile, TotalIntensity from peptidequants");
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
					
					int index=Collections.binarySearch(sourceFiles, sourceFile);
					if (index<0) throw new EncyclopediaException("Unexpected sample: "+sourceFile);
					
					float[] array=intensitiesByPeptideModSeq.get(peptideModSeq);
					if (array==null) {
						array=new float[sourceFiles.size()];
						intensitiesByPeptideModSeq.put(peptideModSeq, array);
					}
					float tic=ticBySourceFileMap.get(sourceFile);
					if (tic>0) {
						array[index]+=totalIntensity;///tic; // sums charge states together
					}
				}
				Logger.logLine("Finished processing "+count+" records, writing report...");
				
				for (Entry<String, float[]> entry : intensitiesByPeptideModSeq.entrySet()) {
					writer.print(entry.getKey());
					float[] array=entry.getValue();
					
					if (false&&array.length>1) {
						float mean=General.mean(array);
						float stdev=General.stdev(array);
						for (int i=0; i<array.length; i++) {
							if (stdev==0) {
								array[i]=0.0f;
							} else {
								array[i]=(array[i]-mean)/stdev;
							}
						}
					}
					
					
					for (float f : array) {
						writer.print("\t"+f);
					}
					writer.println();
				}
				Logger.logLine("Finished writing report!");
				
				rs.close();
			} finally {
				s.close();
				if (writer!=null) writer.close();
			}
		} finally {
			c.close();
		}
	}
}
