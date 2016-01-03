package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class BlibFile extends SQLFile {
	private final File tempFile;
	private File userFile;

	public BlibFile() throws IOException {
		tempFile=File.createTempFile("encyclopedia_", ".blib");
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
	
	public void getStreamEntriesToLibrary(LibraryFile library) throws IOException, SQLException, DataFormatException {
		library.dropIndices();
		
		Connection c=getConnection(userFile);
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s
						.executeQuery("select RefSpectra.precursorMZ, RefSpectra.precursorCharge, RefSpectra.peptideModSeq, RefSpectra.copies, RefSpectra.numPeaks, RefSpectra.retentionTime, RefSpectra.score, RefSpectraPeaks.peakMZ, RefSpectraPeaks.peakIntensity from RefSpectra, RefSpectraPeaks "
								+"where RefSpectra.id == RefSpectraPeaks.RefSpectraID");

				ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
				while (rs.next()) {
					double precursorMZ=rs.getDouble(1);
					byte precursorCharge=(byte)rs.getInt(2);
					String peptideModSeq=rs.getString(3);
					int copies=rs.getInt(4);
					int numPeaks=rs.getInt(5);
					float retentionTime=(float)rs.getDouble(6);
					float score=(float)rs.getDouble(7);
					double[] massArray=decompressDouble(rs.getBytes(8), numPeaks);
					float[] intensityArray=decompressFloat(rs.getBytes(9), numPeaks);

					entries.add(new LibraryEntry(precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray));
					
					if (entries.size()>1000) {
						library.addEntries(entries);
						entries.clear();
						Logger.log(".");
					}
				}

				library.addEntries(entries);
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
		
		library.createIndices();
	}

	private double[] decompressDouble(byte[] bytes, int numPeaks) throws IOException, DataFormatException {
		byte[] decompressed;
		if (bytes.length==numPeaks*8) {
			decompressed=bytes;
		} else {
			decompressed=CompressionUtils.decompress(bytes, numPeaks*8);
		}
		return ByteConverter.toDoubleArray(decompressed, ByteOrder.LITTLE_ENDIAN);
	}

	private float[] decompressFloat(byte[] bytes, int numPeaks) throws IOException, DataFormatException {
		byte[] decompressed;
		if (bytes.length==numPeaks*4) {
			decompressed=bytes;
		} else {
			decompressed=CompressionUtils.decompress(bytes, numPeaks*4);
		}
		return ByteConverter.toFloatArray(decompressed, ByteOrder.LITTLE_ENDIAN);
	}
	
	private void createNewTables() throws IOException, SQLException {
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				s.execute("CREATE TABLE LibInfo(libLSID TEXT, createTime TEXT, numSpecs INTEGER, majorVersion INTEGER, minorVersion INTEGER");
				s.execute("CREATE TABLE Modifications (id INTEGER primary key autoincrement not null,RefSpectraID INTEGER, position INTEGER, mass REAL");
				s.execute("CREATE TABLE RefSpectra (id INTEGER primary key autoincrement not null, peptideSeq VARCHAR(150), precursorMZ REAL, precursorCharge INTEGER, peptideModSeq VARCHAR(200), prevAA CHAR(1), nextAA CHAR(1), copies INTEGER, numPeaks INTEGER, ionMobilityValue REAL, ionMobilityType INTEGER, retentionTime REAL, fileID INTEGER, SpecIDinFile VARCHAR(256), score REAL, scoreType TINYINT");
				s.execute("CREATE TABLE RefSpectraPeaks(RefSpectraID INTEGER, peakMZ BLOB, peakIntensity BLOB");
				s.execute("CREATE TABLE RetentionTimes (RefSpectraID INTEGER, RedundantRefSpectraID INTEGER, SpectrumSourceID INTEGER, ionMobilityValue REAL, ionMobilityType INTEGER, retentionTime REAL, bestSpectrum INTEGER, FOREIGN KEY(RefSpectraID) REFERENCES RefSpectra(id) ");
				s.execute("CREATE TABLE ScoreTypes (id INTEGER PRIMARY KEY, scoreType VARCHAR(128) ");
				s.execute("CREATE TABLE SpectrumSourceFiles (id INTEGER PRIMARY KEY autoincrement not null,fileName VARCHAR(512) ");
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
				s.execute("drop index if exists idxPeptideMod");
				s.execute("drop index if exists idxRefIdPeaks");

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
				s.execute("CREATE INDEX idxPeptide ON RefSpectra (peptideSeq, precursorCharge)");
				s.execute("CREATE INDEX idxPeptideMod ON RefSpectra (peptideModSeq, precursorCharge)");
				s.execute("CREATE INDEX idxRefIdPeaks ON RefSpectraPeaks (RefSpectraID)");
				
				c.commit();
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}
}
