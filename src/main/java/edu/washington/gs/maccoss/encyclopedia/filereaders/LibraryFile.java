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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class LibraryFile extends SQLFile implements LibraryInterface {
	public static final String ELIB=".elib";
	public static final String VERSION_STRING="version";
	public static final String MOST_RECENT_VERSION="0.1";
	
	private File userFile=null;
	private final File tempFile;

	public LibraryFile() throws IOException {
		tempFile=File.createTempFile("encyclopedia_", ELIB);
		tempFile.deleteOnExit();
	}
	
	public String getName() {
		return userFile==null?tempFile.getName():userFile.getName();
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

	public void saveAsFile(File userFile) throws IOException, SQLException {
		this.userFile=userFile;
		saveFile();
	}

	public void saveFile() throws IOException, SQLException {
		if (userFile!=null) {
			setFileVersion();
			
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

	public void setFileName(String fileName, String sourceName, String fileLocation) throws IOException, SQLException {
		HashMap<String, String> map=new HashMap<String, String>();
		map.put("filename", fileName);
		map.put("sourcename", sourceName);
		map.put("filelocation", fileLocation);
		addMetadata(map);
	}
	
	public void setFileVersion() throws IOException, SQLException {
		HashMap<String, String> map=new HashMap<String, String>();
		map.put(VERSION_STRING, MOST_RECENT_VERSION);
		addMetadata(map);
	}

	public void addMetadata(Map<String, String> data) throws IOException, SQLException {
		Connection c=getConnection(tempFile);
		try {
			PreparedStatement prep=c.prepareStatement("insert into metadata (Key, Value) VALUES (?,?)");
			try {
				for (Entry<String, String> entry : data.entrySet()) {
					prep.setString(1, entry.getKey());
					prep.setString(2, entry.getValue());
					prep.addBatch();
				}
				prep.executeBatch();
				prep.close();
				c.commit();
			} finally {
				prep.close();
			}
		} finally {
			c.close();
		}
	}
	
	public HashMap<String, String> getMetadata() throws IOException, SQLException {
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select Key, Value from metadata");

				HashMap<String, String> map=new HashMap<String, String>();
				while (rs.next()) {
					String key=rs.getString(1);
					String value=rs.getString(2);
					map.put(key, value);
				}

				return map;
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	public void addIntegratedEntries(ArrayList<IntegratedLibraryEntry> entries, float minimumCorrelation) throws IOException, SQLException {
		// first add normal data
		ArrayList<LibraryEntry> recasted=new ArrayList<LibraryEntry>();
		for (IntegratedLibraryEntry entry : entries) {
			recasted.add(entry);
		}
		addEntries(recasted);
		
		// then add integrated areas
		Connection c=getConnection(tempFile);
		try {
			PreparedStatement peptidePrep=c.prepareStatement("INSERT INTO peptidequants (PrecursorCharge, PeptideModSeq, SourceFile, RTInSecondsStart, RTInSecondsStop, TotalIntensity, NumberOfQuantIons, MedianChromatogramEncodedLength, MedianChromatogramArray) VALUES (?,?,?,?,?,?,?,?,?)");
			PreparedStatement fragmentPrep=c.prepareStatement("INSERT INTO fragmentquants (PrecursorCharge, PeptideModSeq, SourceFile, IonType, FragmentMass, Correlation, Intensity) VALUES (?,?,?,?,?,?,?)");
			try {
				for (IntegratedLibraryEntry entry : entries) {
					TransitionRefinementData data=entry.getRefinementData();
					
					peptidePrep.setInt(1, entry.getPrecursorCharge());
					peptidePrep.setString(2, entry.getPeptideModSeq());
					peptidePrep.setString(3, entry.getSource());
					peptidePrep.setFloat(4, data.getRange().getStart());
					peptidePrep.setFloat(5, data.getRange().getStop());
					peptidePrep.setFloat(6, data.getTotalIntensity(minimumCorrelation));
					peptidePrep.setInt(7, data.getTotalQuantIons(minimumCorrelation));
					byte[] intensityByteArray=ByteConverter.toByteArray(data.getMedianChromatogram());
					peptidePrep.setInt(8, intensityByteArray.length);
					peptidePrep.setBytes(9, CompressionUtils.compress(intensityByteArray));
					peptidePrep.addBatch();

					float[] correlationArray=data.getCorrelationArray();
					for (int i=0; i<correlationArray.length; i++) {
						if (correlationArray[i]>=minimumCorrelation) {
							fragmentPrep.setInt(1, entry.getPrecursorCharge());
							fragmentPrep.setString(2, entry.getPeptideModSeq());
							fragmentPrep.setString(3, entry.getSource());
							fragmentPrep.setString(4, "BY");
							fragmentPrep.setDouble(5, data.getFragmentMassArray()[i]);
							fragmentPrep.setFloat(6, correlationArray[i]);
							fragmentPrep.setFloat(7, data.getIntegrationArray()[i]);
							fragmentPrep.addBatch();
						}
					}
					
				}
				peptidePrep.executeBatch();
				fragmentPrep.executeBatch();
		
				c.commit();
			} finally {
				peptidePrep.close();
				fragmentPrep.close();
			}
		} finally {
			c.close();
		}
	}
	
	public void addEntries(ArrayList<LibraryEntry> entries) throws IOException, SQLException {
		Connection c=getConnection(tempFile);
		try {
			PreparedStatement prep=c.prepareStatement("INSERT INTO entries (PrecursorMZ, PrecursorCharge, PeptideModSeq, PeptideSeq, Copies, RTInSeconds, Score, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray, SourceFile) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)");
			PreparedStatement proteinPrep=c.prepareStatement("INSERT OR IGNORE INTO proteins (PeptideSeq, ProteinAccessions) VALUES (?,?)");
			try {
				for (LibraryEntry entry : entries) {
					if (entry.getAccessions().size()==0) continue;
					
					String pepSeq=entry.getPeptideSeq();
					prep.setDouble(1, entry.getPrecursorMZ());
					prep.setInt(2, entry.getPrecursorCharge());
					prep.setString(3, entry.getPeptideModSeq());
					prep.setString(4, pepSeq);
					prep.setInt(5, entry.getCopies());
					prep.setFloat(6, entry.getRetentionTime());
					prep.setFloat(7, entry.getScore());
					byte[] massByteArray=ByteConverter.toByteArray(entry.getMassArray());
					prep.setInt(8, massByteArray.length);
					prep.setBytes(9, CompressionUtils.compress(massByteArray));
					byte[] intensityByteArray=ByteConverter.toByteArray(entry.getIntensityArray());
					prep.setInt(10, intensityByteArray.length);
					prep.setBytes(11, CompressionUtils.compress(intensityByteArray));
					prep.setString(12, entry.getSource());
					prep.addBatch();
					
					proteinPrep.setString(1, pepSeq);
					proteinPrep.setString(2, PSMData.accessionsToString(entry.getAccessions()));
					proteinPrep.addBatch();
				}
				prep.executeBatch();
				proteinPrep.executeBatch();
				
				c.commit();
			} finally {
				prep.close();
				proteinPrep.close();
			}
		} finally {
			c.close();
		}
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface#getEntries(edu.washington.gs.maccoss.encyclopedia.datastructures.Range)
	 */
	@Override
	public ArrayList<LibraryEntry> getEntries(String peptideModSeq, byte charge, boolean sqrt) throws IOException, SQLException, DataFormatException {
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select e.PrecursorMZ, e.PrecursorCharge, e.PeptideModSeq, e.Copies, e.RTInSeconds, e.Score, e.MassEncodedLength, e.MassArray, e.IntensityEncodedLength, e.IntensityArray, p.ProteinAccessions, e.SourceFile from entries e, proteins p"
						+ " where e.PeptideSeq=p.PeptideSeq and e.PeptideModSeq = \""+peptideModSeq+"\" and e.PrecursorCharge = "+charge);

				ArrayList<LibraryEntry> entry=new ArrayList<LibraryEntry>();
				while (rs.next()) {

					double precursorMZ=rs.getDouble(1);
					byte precursorCharge=(byte)rs.getInt(2);
					peptideModSeq=rs.getString(3);
					int copies=rs.getInt(4);
					float retentionTime=rs.getFloat(5); 
					float score=rs.getFloat(6);
					int massEncodedLength=rs.getInt(7);
					double[] massArray=ByteConverter.toDoubleArray(CompressionUtils.decompress(rs.getBytes(8), massEncodedLength));
					int intensityEncodedLength=rs.getInt(9);
					float[] intensityArray=ByteConverter.toFloatArray(CompressionUtils.decompress(rs.getBytes(10), intensityEncodedLength));
					if (sqrt) {
						intensityArray=General.protectedSqrt(intensityArray);
					}
					HashSet<String> accessions=PSMData.stringToAccessions(rs.getString(11));
					String sourceFile=rs.getString(12);
					entry.add(new LibraryEntry(sourceFile, accessions, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray));
				}

				return entry;
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface#getEntries(edu.washington.gs.maccoss.encyclopedia.datastructures.Range)
	 */
	@Override
	public Range getMinMaxMZ() throws IOException, SQLException {
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select min(PrecursorMZ), max(PrecursorMZ) from entries");

				while (rs.next()) {
					double min=rs.getDouble(1);
					double max=rs.getDouble(2);
					return new Range((float)min, (float)max);
				}

				return new Range(0, 0);
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface#getEntries(edu.washington.gs.maccoss.encyclopedia.datastructures.Range)
	 */
	@Override
	public ArrayList<LibraryEntry> getEntries(Range precursorMz, boolean sqrt) throws IOException, SQLException, DataFormatException {
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select e.PrecursorMZ, e.PrecursorCharge, e.PeptideModSeq, e.Copies, e.RTInSeconds, e.Score, e.MassEncodedLength, e.MassArray, e.IntensityEncodedLength, e.IntensityArray, p.ProteinAccessions, e.SourceFile from entries e, proteins p"
						+ " where e.PeptideSeq=p.PeptideSeq and e.PrecursorMz between "+precursorMz.getStart()+" and "+precursorMz.getStop());

				ArrayList<LibraryEntry> entry=new ArrayList<LibraryEntry>();
				while (rs.next()) {

					double precursorMZ=rs.getDouble(1);
					byte precursorCharge=(byte)rs.getInt(2);
					String peptideModSeq=rs.getString(3);
					int copies=rs.getInt(4);
					float retentionTime=rs.getFloat(5); 
					float score=rs.getFloat(6);
					int massEncodedLength=rs.getInt(7);
					double[] massArray=ByteConverter.toDoubleArray(CompressionUtils.decompress(rs.getBytes(8), massEncodedLength));
					int intensityEncodedLength=rs.getInt(9);
					float[] intensityArray=ByteConverter.toFloatArray(CompressionUtils.decompress(rs.getBytes(10), intensityEncodedLength));
					if (sqrt) {
						intensityArray=General.protectedSqrt(intensityArray);
					}
					HashSet<String> accessions=PSMData.stringToAccessions(rs.getString(11));
					String sourceFile=rs.getString(12);
					entry.add(new LibraryEntry(sourceFile, accessions, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray));
				}

				return entry;
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}
	
	public ArrayList<LibraryEntry> getAllEntries(boolean sqrt) throws IOException, SQLException, DataFormatException {
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select e.PrecursorMZ, e.PrecursorCharge, e.PeptideModSeq, e.Copies, e.RTInSeconds, e.Score, e.MassEncodedLength, e.MassArray, e.IntensityEncodedLength, e.IntensityArray, p.ProteinAccessions, e.SourceFile from entries e, proteins p"
						+ " where e.PeptideSeq=p.PeptideSeq");

				ArrayList<LibraryEntry> entry=new ArrayList<LibraryEntry>();
				while (rs.next()) {

					double precursorMZ=rs.getDouble(1);
					byte precursorCharge=(byte)rs.getInt(2);
					String peptideModSeq=rs.getString(3);
					int copies=rs.getInt(4);
					float retentionTime=rs.getFloat(5); 
					float score=rs.getFloat(6);
					int massEncodedLength=rs.getInt(7);
					double[] massArray=ByteConverter.toDoubleArray(CompressionUtils.decompress(rs.getBytes(8), massEncodedLength));
					int intensityEncodedLength=rs.getInt(9);
					float[] intensityArray=ByteConverter.toFloatArray(CompressionUtils.decompress(rs.getBytes(10), intensityEncodedLength));
					if (sqrt) {
						intensityArray=General.protectedSqrt(intensityArray);
					}
					HashSet<String> accessions=PSMData.stringToAccessions(rs.getString(11));
					String sourceFile=rs.getString(12);
					entry.add(new LibraryEntry(sourceFile, accessions, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray));
				}

				return entry;
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	private void createNewTables() throws IOException, SQLException {
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				s.execute("create table if not exists metadata ( Key string not null, Value string not null, primary key (Key) )");
				s.execute("create table if not exists entries ( PrecursorMz double not null, PrecursorCharge int not null, PeptideModSeq string not null, PeptideSeq string not null, Copies int not null, RTInSeconds double not null, Score double not null, MassEncodedLength int not null, MassArray blob not null, IntensityEncodedLength int not null, IntensityArray blob not null, SourceFile string not null )");
				s.execute("create table if not exists proteins ( PeptideSeq string not null, ProteinAccessions string not null, primary key (PeptideSeq) )");
				s.execute("create table if not exists peptidequants ( PrecursorCharge int not null, PeptideModSeq string not null, SourceFile string not null, RTInSecondsStart double not null, RTInSecondsStop double not null, TotalIntensity double not null, NumberOfQuantIons int not null, MedianChromatogramEncodedLength int not null, MedianChromatogramArray blob not null )");
				s.execute("create table if not exists fragmentquants ( PrecursorCharge int not null, PeptideModSeq string not null, SourceFile string not null, IonType string not null, FragmentMass double not null, Correlation double not null, Intensity double not null )");
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
				
				s.execute("drop index if exists \"PeptideModSeq_Entries_index\"");
				s.execute("drop index if exists \"PeptideSeq_Entries_index\"");
				s.execute("drop index if exists \"PrecursorMz_Entries_index\"");
				s.execute("drop index if exists \"ProteinAccessions_Proteins_index\"");
				s.execute("drop index if exists \"PeptideModSeq_Peptides_index\"");
				s.execute("drop index if exists \"PeptideModSeq_Fragments_index\"");

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
				
				s.execute("create index if not exists \"PeptideModSeq_Entries_index\" on \"entries\" (\"PeptideModSeq\" ASC)");
				s.execute("create index if not exists \"PeptideSeq_Entries_index\" on \"entries\" (\"PeptideSeq\" ASC)");
				s.execute("create index if not exists \"PrecursorMz_Entries_index\" on \"entries\" (\"PrecursorMz\" ASC)");
				s.execute("create index if not exists \"ProteinAccessions_Proteins_index\" on \"proteins\" (\"ProteinAccessions\" ASC)");
				s.execute("create index if not exists \"PeptideModSeq_Peptides_index\" on \"peptidequants\" (\"PeptideModSeq\" ASC)");
				s.execute("create index if not exists \"PeptideModSeq_Fragments_index\" on \"fragmentquants\" (\"PeptideModSeq\" ASC)");

				c.commit();
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	public void close() {
		if (!tempFile.delete()) {
			Logger.errorLine("Error deleting temp file!");
		}
	}

}
