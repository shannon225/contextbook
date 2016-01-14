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
import java.util.Map;
import java.util.Map.Entry;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class LibraryFile extends SQLFile {
	private File userFile=null;
	private final File tempFile;

	public LibraryFile() throws IOException {
		tempFile=File.createTempFile("encyclopedia_", ".elib");
		tempFile.deleteOnExit();
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

	public void addPrecursor(ArrayList<PrecursorScan> precursors) throws IOException, SQLException {
		Connection c=getConnection(tempFile);
		try {
			PreparedStatement prep=c.prepareStatement("insert into precursor (SpectrumName, SpectrumIndex, ScanStartTime, PeakCount, MassArray, IntensityArray) VALUES (?,?,?,?,?,?)");
			try {
				for (PrecursorScan precursor : precursors) {
					prep.setString(1, precursor.getSpectrumName());
					prep.setInt(2, precursor.getSpectrumIndex());
					prep.setFloat(3, precursor.getScanStartTime());
					prep.setFloat(4, precursor.getMassArray().length);
					prep.setBytes(5, CompressionUtils.compress(ByteConverter.toByteArray(precursor.getMassArray())));
					prep.setBytes(6, CompressionUtils.compress(ByteConverter.toByteArray(precursor.getIntensityArray())));
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

	public void addEntries(ArrayList<LibraryEntry> entries) throws IOException, SQLException {
		Connection c=getConnection(tempFile);
		try {
			PreparedStatement prep=c
					.prepareStatement("insert into entries (PrecursorMZ, PrecursorCharge, PeptideModSeq, Copies, RetentionTime, Score, PeakCount, MassArray, IntensityArray) VALUES (?,?,?,?,?,?,?,?,?)");
			try {
				for (LibraryEntry entry : entries) {
					prep.setDouble(1, entry.getPrecursorMZ());
					prep.setInt(2, entry.getPrecursorCharge());
					prep.setString(3, entry.getPeptideModSeq());
					prep.setInt(4, entry.getCopies());
					prep.setFloat(5, entry.getRetentionTime());
					prep.setFloat(6, entry.getScore());
					prep.setFloat(7, entry.getMassArray().length);
					prep.setBytes(8, CompressionUtils.compress(ByteConverter.toByteArray(entry.getMassArray())));
					prep.setBytes(9, CompressionUtils.compress(ByteConverter.toByteArray(entry.getIntensityArray())));
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

	public ArrayList<LibraryEntry> getEntries(Range precursorMz) throws IOException, SQLException, DataFormatException {
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select PrecursorMZ, PrecursorCharge, PeptideModSeq, Copies, RetentionTime, Score, PeakCount, MassArray, IntensityArray from entries " +
						"where PrecursorMz between "+precursorMz.getStart()+" and "+precursorMz.getStop());

				ArrayList<LibraryEntry> entry=new ArrayList<LibraryEntry>();
				while (rs.next()) {

					double precursorMZ=rs.getDouble(1);
					byte precursorCharge=(byte)rs.getInt(2);
					String peptideModSeq=rs.getString(3);
					int copies=rs.getInt(4);
					float retentionTime=rs.getFloat(5); 
					float score=rs.getFloat(6); 
					int peakCount=rs.getInt(7);
					double[] massArray=ByteConverter.toDoubleArray(CompressionUtils.decompress(rs.getBytes(8), peakCount));
					float[] intensityArray=ByteConverter.toFloatArray(CompressionUtils.decompress(rs.getBytes(9), peakCount));
					entry.add(new LibraryEntry(precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray));
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
				s.execute("create table if not exists entries (PrecursorMz double not null, PrecursorCharge int not null, PeptideModSeq string not null, Copies int not null, RetentionTime double not null, Score double not null, PeakCount, MassArray blob not null, IntensityArray blob not null)");

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
				
				s.execute("drop index if exists \"PeptideModSeq_index\"");
				s.execute("drop index if exists \"PrecursorMz_index\"");

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
				
				s.execute("create index if not exists \"PeptideModSeq_index\" on \"entries\" (\"PeptideModSeq\" ASC)");
				s.execute("create index if not exists \"PrecursorMz_index\" on \"entries\" (\"PrecursorMz\" ASC)");

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
