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

import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class StripeFile extends SQLFile implements StripeFileInterface {
	public static final String DIA_EXTENSION=".dia";
	
	private File userFile=null;
	private final File tempFile;
	
	private final HashMap<Range, Float> ranges=new HashMap<Range, Float>();

	public StripeFile() throws IOException {
		tempFile=File.createTempFile("encyclopedia_", DIA_EXTENSION);
		tempFile.deleteOnExit();
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface#getRanges()
	 */
	@Override
	@SuppressWarnings("unchecked")
	public HashMap<Range, Float> getRanges() {
		return (HashMap<Range, Float>)ranges.clone();
	}
	
	public void setRanges(HashMap<Range, Float> ranges) {
		this.ranges.clear();
		this.ranges.putAll(ranges);
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface#openFile(java.io.File)
	 */
	@Override
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
				ResultSet rs=s.executeQuery("select Start, Stop, DutyCycle from Ranges");

				while (rs.next()) {
					float start=rs.getFloat(1);
					float stop=rs.getFloat(2);
					float dutyCycle=rs.getFloat(3);
					ranges.put(new Range(start, stop), dutyCycle);
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
			PreparedStatement prep=c.prepareStatement("insert into ranges (Start, Stop, DutyCycle) VALUES (?,?,?)");
			try {
				for (Entry<Range, Float> entry : ranges.entrySet()) {
					Range range=entry.getKey();
					prep.setFloat(1, range.getStart());
					prep.setFloat(2, range.getStop());
					prep.setFloat(3, entry.getValue());
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

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface#openFile()
	 */
	@Override
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
			PreparedStatement prep=c.prepareStatement("insert into precursor (SpectrumName, SpectrumIndex, ScanStartTime, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray) VALUES (?,?,?,?,?,?,?)");
			try {
				for (PrecursorScan precursor : precursors) {
					prep.setString(1, precursor.getSpectrumName());
					prep.setInt(2, precursor.getSpectrumIndex());
					prep.setFloat(3, precursor.getScanStartTime());
					byte[] massByteArray=ByteConverter.toByteArray(precursor.getMassArray());
					prep.setInt(4, massByteArray.length);
					prep.setBytes(5, CompressionUtils.compress(massByteArray));
					byte[] intensityByteArray=ByteConverter.toByteArray(precursor.getIntensityArray());
					prep.setInt(6, intensityByteArray.length);
					prep.setBytes(7, CompressionUtils.compress(intensityByteArray));
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

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface#getPrecursors(float, float)
	 */
	@Override
	public ArrayList<PrecursorScan> getPrecursors(float minRT, float maxRT) throws IOException, SQLException,DataFormatException {
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select SpectrumName, SpectrumIndex, ScanStartTime, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray from precursor "
						+"where ScanStartTime between "+minRT+" and "+maxRT);

				ArrayList<PrecursorScan> precursors=new ArrayList<PrecursorScan>();
				while (rs.next()) {
					String spectrumName=rs.getString(1);
					int spectrumIndex=rs.getInt(2);
					float scanStartTime=rs.getFloat(3);
					int massEncodedLength=rs.getInt(4);
					double[] massArray=ByteConverter.toDoubleArray(CompressionUtils.decompress(rs.getBytes(5), massEncodedLength));
					int intensityEncodedLength=rs.getInt(6);
					float[] intensityArray=ByteConverter.toFloatArray(CompressionUtils.decompress(rs.getBytes(7), intensityEncodedLength));
					precursors.add(new PrecursorScan(spectrumName, spectrumIndex, scanStartTime, massArray, intensityArray));
				}

				return precursors;
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	public void addStripe(ArrayList<Stripe> stripes) throws IOException, SQLException {
		Connection c=getConnection(tempFile);
		try {
			PreparedStatement prep=c
					.prepareStatement("insert into spectra (SpectrumName, PrecursorName, SpectrumIndex, ScanStartTime, IsolationWindowLower, IsolationWindowUpper, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray) VALUES (?,?,?,?,?,?,?,?,?,?)");
			try {
				for (Stripe stripe : stripes) {
					prep.setString(1, stripe.getSpectrumName());
					prep.setString(2, stripe.getPrecursorName());
					prep.setInt(3, stripe.getSpectrumIndex());
					prep.setFloat(4, stripe.getScanStartTime());
					prep.setFloat(5, stripe.getIsolationWindowLower());
					prep.setFloat(6, stripe.getIsolationWindowUpper());
					byte[] massByteArray=ByteConverter.toByteArray(stripe.getMassArray());
					prep.setInt(7, massByteArray.length);
					prep.setBytes(8, CompressionUtils.compress(massByteArray));
					byte[] intensityByteArray=ByteConverter.toByteArray(stripe.getIntensityArray());
					prep.setInt(9, intensityByteArray.length);
					prep.setBytes(10, CompressionUtils.compress(intensityByteArray));
					
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

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface#getStripes(double, float, float, boolean)
	 */
	@Override
	public ArrayList<Stripe> getStripes(double targetMz, float minRT, float maxRT, boolean sqrt) throws IOException, SQLException,DataFormatException {
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select SpectrumName, PrecursorName, SpectrumIndex, ScanStartTime, IsolationWindowLower, IsolationWindowUpper, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray from spectra "
						+"where IsolationWindowLower <= "+targetMz+" and IsolationWindowUpper >= "+targetMz+" and ScanStartTime between "+minRT+" and "+maxRT);

				ArrayList<Stripe> stripes=new ArrayList<Stripe>();
				while (rs.next()) {
					String spectrumName=rs.getString(1);
					String precursorName=rs.getString(2);
					int spectrumIndex=rs.getInt(3);
					float scanStartTime=rs.getFloat(4);
					float isolationWindowLower=rs.getFloat(5);
					float isolationWindowUpper=rs.getFloat(6);
					int massEncodedLength=rs.getInt(7);
					double[] massArray=ByteConverter.toDoubleArray(CompressionUtils.decompress(rs.getBytes(8), massEncodedLength));
					int intensityEncodedLength=rs.getInt(9);
					float[] intensityArray=ByteConverter.toFloatArray(CompressionUtils.decompress(rs.getBytes(10), intensityEncodedLength));
					if (sqrt) {
						intensityArray=General.protectedSqrt(intensityArray);
					}
					stripes.add(new Stripe(spectrumName, precursorName, spectrumIndex, scanStartTime, isolationWindowLower, isolationWindowUpper, massArray, intensityArray));
				}

				return stripes;
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
				s.execute("create table if not exists ranges ( Start float not null, Stop float not null, DutyCycle float not null )");
				s.execute("create table if not exists spectra ( SpectrumName string not null, PrecursorName string not null, SpectrumIndex int not null, ScanStartTime float not null, IsolationWindowLower float not null, IsolationWindowUpper float not null, MassEncodedLength int not null, MassArray blob not null, IntensityEncodedLength int not null, IntensityArray blob not null, primary key (SpectrumIndex) )");
				s.execute("create table if not exists precursor ( SpectrumName string not null, SpectrumIndex int not null, ScanStartTime float not null, MassEncodedLength int not null, MassArray blob not null, IntensityEncodedLength int not null, IntensityArray blob not null, primary key (SpectrumIndex) )");

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

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface#close()
	 */
	@Override
	public void close() {
		if (!tempFile.delete()) {
			Logger.errorLine("Error deleting temp file!");
		}
	}
}
