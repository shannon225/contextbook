package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.nio.ByteOrder;
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
	private final File userFile;

	public BlibFile(File userFile) {
		this.userFile=userFile;
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
					int precursorCharge=rs.getInt(2);
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
}
