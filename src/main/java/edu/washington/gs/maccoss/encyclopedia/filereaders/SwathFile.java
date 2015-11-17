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

import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Swath;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class SwathFile extends SQLFile {
	private File userFile=null;
	private final File tempFile;
	
	private final HashSet<Range> ranges=new HashSet<Range>();

	public SwathFile() throws IOException {
		tempFile=File.createTempFile("encyclopedia_", ".dia");
		tempFile.deleteOnExit();
	}
	
	@SuppressWarnings("unchecked")
	public HashSet<Range> getRanges() {
		return (HashSet<Range>)ranges.clone();
	}

	public void openFile(File userFile) throws IOException, SQLException {
		this.userFile=userFile;
		openFile();
		loadRanges();
	}
	
	public void loadRanges() throws IOException, SQLException {
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select Start, Stop from Ranges");

				while (rs.next()) {
					float start=rs.getFloat(1);
					float stop=rs.getFloat(2);
					ranges.add(new Range(start, stop));
				}
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}
	
	public void writeRanges() throws IOException, SQLException {
		Connection c=getConnection(tempFile);
		try {
			PreparedStatement prep=c.prepareStatement("insert into ranges (Start, Stop) VALUES (?,?)");
			try {
				for (Range entry : ranges) {
					prep.setFloat(1, entry.getStart());
					prep.setFloat(2, entry.getStop());
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
		writeRanges();
		
		if (userFile!=null) {
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
					prep.setInt(4, precursor.getMassArray().length);
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

	public void addSwath(ArrayList<Swath> swaths) throws IOException, SQLException {
		Connection c=getConnection(tempFile);
		try {
			PreparedStatement prep=c
					.prepareStatement("insert into spectra (SpectrumName, PrecursorName, SpectrumIndex, ScanStartTime, IsolationWindowLower, IsolationWindowUpper, PeakCount, MassArray, IntensityArray) VALUES (?,?,?,?,?,?,?,?,?)");
			try {
				for (Swath swath : swaths) {
					ranges.add(swath.getRange());
					
					prep.setString(1, swath.getSpectrumName());
					prep.setString(2, swath.getPrecursorName());
					prep.setInt(3, swath.getSpectrumIndex());
					prep.setFloat(4, swath.getScanStartTime());
					prep.setFloat(5, swath.getIsolationWindowLower());
					prep.setFloat(6, swath.getIsolationWindowUpper());
					prep.setInt(7, swath.getMassArray().length);
					prep.setBytes(8, CompressionUtils.compress(ByteConverter.toByteArray(swath.getMassArray())));
					prep.setBytes(9, CompressionUtils.compress(ByteConverter.toByteArray(swath.getIntensityArray())));
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

	public ArrayList<Swath> getSwaths(double targetMz, float minRT, float maxRT, boolean sqrt) throws IOException, SQLException,DataFormatException {
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select SpectrumName, PrecursorName, SpectrumIndex, ScanStartTime, IsolationWindowLower, IsolationWindowUpper, PeakCount, MassArray, IntensityArray from spectra "
						+"where IsolationWindowLower <= "+targetMz+" and IsolationWindowUpper >= "+targetMz+" and ScanStartTime between "+minRT+" and "+maxRT);

				ArrayList<Swath> swaths=new ArrayList<Swath>();
				while (rs.next()) {
					String spectrumName=rs.getString(1);
					String precursorName=rs.getString(2);
					int spectrumIndex=rs.getInt(3);
					float scanStartTime=rs.getFloat(4);
					float isolationWindowLower=rs.getFloat(5);
					float isolationWindowUpper=rs.getFloat(6);
					int peakCount=rs.getInt(7);
					double[] massArray=ByteConverter.toDoubleArray(CompressionUtils.decompress(rs.getBytes(8), peakCount));
					float[] intensityArray=ByteConverter.toFloatArray(CompressionUtils.decompress(rs.getBytes(9), peakCount));
					if (sqrt) {
						intensityArray=General.protectedSqrt(intensityArray);
					}
					swaths.add(new Swath(spectrumName, precursorName, spectrumIndex, scanStartTime, isolationWindowLower, isolationWindowUpper, massArray, intensityArray));
				}

				return swaths;
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
				s.execute("create table if not exists ranges ( Start float not null, Stop float not null )");
				s.execute("create table if not exists spectra ( SpectrumName string not null, PrecursorName string not null, SpectrumIndex int not null, ScanStartTime float not null, IsolationWindowLower float not null, IsolationWindowUpper float not null, PeakCount int not null, MassArray blob not null, IntensityArray blob not null, primary key (SpectrumIndex) )");
				s.execute("create table if not exists precursor ( SpectrumName string not null, SpectrumIndex int not null, ScanStartTime float not null, PeakCount int not null, MassArray blob not null, IntensityArray blob not null, primary key (SpectrumIndex) )");

				s.execute("create index if not exists \"spectra_index_isolation_window_lower\" on \"spectra\" (\"IsolationWindowLower\" ASC)");
				s.execute("create index if not exists \"spectra_index_isolation_window_upper\" on \"spectra\" (\"IsolationWindowUpper\" ASC)");
				s.execute("create index if not exists \"spectra_index_scan_start_time\" on \"spectra\" (\"ScanStartTime\" ASC)");

				s.execute("create index if not exists \"precursor_index_scan_start_time\" on \"precursor\" (\"ScanStartTime\" ASC)");

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
