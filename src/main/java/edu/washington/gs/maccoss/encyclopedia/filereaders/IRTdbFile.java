package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ChromatogramLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class IRTdbFile extends SQLFile {
	public static final String IRTDB=".irtdb";
	private final File tempFile;
	private File userFile;

	public IRTdbFile() throws IOException {
		tempFile=File.createTempFile("encyclopedia_", IRTDB);
		tempFile.deleteOnExit();
	}
	public void addLibrary(ArrayList<LibraryEntry> entries) throws IOException, SQLException, DataFormatException {

		Connection c=getConnection(tempFile);
		try {
			Statement normalStatement=c.createStatement();
			try {
				normalStatement.executeUpdate("insert into VersionInfo (SchemaVersion) VALUES (1);");
				//ResultSet rs=normalStatement.executeQuery("select seq from sqlite_sequence");
				//int idCounter=rs.getInt(1);
				
				//normalStatement.executeUpdate("insert into sqlite_sequence (name, seq) VALUES (IrtLibrary,\""+(idCounter+entries.size())+"\");");
				
				PreparedStatement prepStatement=c.prepareStatement("insert into IrtLibrary (Id, PeptideModSeq, Irt, Standard, TimeSource) VALUES (?,?,?,?,?)");
				
				// if ChromatogramLibraryEntry, keep the top 10 best entries (based on number of high quality fragment ions)
				// higher scores are better
				HashMap<String, ScoredObject<PeptidePrecursor>> bestEntries=new HashMap<String, ScoredObject<PeptidePrecursor>>();
				for (LibraryEntry entry : entries) {
					if (entry instanceof ChromatogramLibraryEntry) {
						float[] intensityArray = entry.getIntensityArray();
						float intensityThreshold=General.max(intensityArray)*0.05f; // 5% intensity minimum
						float[] correlationArray=entry.getCorrelationArray();
						float score=0.0f;
						int count=0;
						for (int i = 0; i < correlationArray.length; i++) {
							if (correlationArray[i]>TransitionRefiner.quantitativeCorrelationThreshold&&intensityArray[i]>intensityThreshold) {
								count++;
								score+=Log.log10(intensityArray[i]);
								//score+=correlationArray[i];
							}
						}
						// strongly drop peptides with fewer than 6 good transitions (but keep them in if there aren't many options)
						if (count<6) {
							score=score/1000.0f;
						}
						
						ScoredObject<PeptidePrecursor> scored=new ScoredObject<PeptidePrecursor>(score, entry);
						
						ScoredObject<PeptidePrecursor> previous=bestEntries.get(entry.getPeptideModSeq());
						if (previous==null||previous.getScore()<score) {
							// add if new or replace if better
							bestEntries.put(entry.getPeptideModSeq(), scored);
						}
					}
				}
				ArrayList<ScoredObject<PeptidePrecursor>> bestEntriesList=new ArrayList<>(bestEntries.values());
				Collections.sort(bestEntriesList);
				HashSet<String> bestPeptideModSeqs=new HashSet<>();
				int start=Math.max(0, bestEntriesList.size()-10);
				for (int i = start; i < bestEntriesList.size(); i++) {
					bestPeptideModSeqs.add(bestEntriesList.get(i).y.getPeptideModSeq());
					Logger.logLine("Adding "+bestEntriesList.get(i).y.getPeptideModSeq()+" to anchor iRT peptides");
				}

				int idCounter=0;
				for (LibraryEntry entry : entries) {
					// if any possible standards then use the standards, otherwise set everything to be a standard
					boolean isStandard=true;
					if (bestPeptideModSeqs.size()>0) {
						isStandard=bestPeptideModSeqs.contains(entry.getPeptideModSeq());
					}
					idCounter++;
					prepStatement.setInt(1, idCounter); // id
					prepStatement.setString(2, entry.getPeptideSeq()); // pepSeq
					prepStatement.setDouble(3, entry.getScanStartTime()/60f); // rt in min
					prepStatement.setBoolean(4, isStandard);
					prepStatement.setInt(5, isStandard?1:0);
					prepStatement.addBatch();
				}
				prepStatement.executeBatch();
				prepStatement.close();
				
				c.commit();
			} finally {
				normalStatement.close();
			}
		} finally {
			c.close();
		}
	}

	public void close() {
		if (tempFile.exists()&&!tempFile.delete()) {
			Logger.errorLine("Error deleting temp BLIB file!");
		}
	}

	public void openFile(File userFile) throws IOException, SQLException {
		this.userFile=userFile;
		openFile();
	}

	public void openFile() throws IOException, SQLException {
		if (userFile!=null) {
			Files.copy(userFile.toPath(), tempFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
		}
		createNewTables();
	}
	
	public void setUserFile(File userFile) {
		this.userFile=userFile;
	}

	public File getUserFile() {
		return userFile;
	}

	
	public TObjectFloatHashMap<String> getIRTs() throws IOException, SQLException, DataFormatException {
		TObjectFloatHashMap<String> map=new TObjectFloatHashMap<String>();
		
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select PeptideModSeq, Irt from IrtLibrary");
				while (rs.next()) {
					String peptideModSeq=rs.getString(1);
					float retentionTime=(float)rs.getDouble(2);
					map.put(peptideModSeq, retentionTime);
				}
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
		
		return map;
	}

	
	private void createNewTables() throws IOException, SQLException {
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				s.execute("CREATE TABLE IrtLibrary (Id  integer primary key autoincrement, PeptideModSeq TEXT not null, Irt DOUBLE, Standard BOOL, TimeSource INT)");
				s.execute("CREATE TABLE VersionInfo (SchemaVersion INT not null, primary key (SchemaVersion))");
				//s.execute("CREATE TABLE sqlite_sequence(name,seq)");
				c.commit();
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}
	
	public void dropIndices() throws IOException, SQLException {
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				
				s.execute("drop index if exists idxPeptide");

				c.commit();
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	public void createIndices() throws IOException, SQLException {
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				s.execute("CREATE INDEX idxPeptide on IrtLibrary (PeptideModSeq)");
				
				c.commit();
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	public void saveAsFile(File userFile) throws IOException, SQLException {
		this.userFile=userFile;
		saveFile();
	}

	public void saveFile() throws IOException, SQLException {
		if (userFile!=null) {
			Connection c=getConnection(tempFile);

			try {
				Statement s=c.createStatement();
				try {
					s.execute("END");
					s.execute("VACUUM");
					s.execute("BEGIN");
					
					c.commit();
				} finally {
					s.close();
				}
			} finally {
				c.close();
			}
			
			Files.copy(tempFile.toPath(), userFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
		}
	}
}
