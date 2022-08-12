package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class DIAtoMSMSProducer implements MSMSProducer {
	public static final int MAX_PRECURSORS_PER_BLOCK = 1000;
	public static final int MAX_STRIPES_PER_SCAN = 1000;

	private final File diaFile;
	private LibraryFile library; 
	private final BlockingQueue<MSMSBlock> outputBlockQueue;
	private final SearchParameters parameters;
	final private int fraction;

	private Throwable error;

	public DIAtoMSMSProducer(File diaFile, int fraction, BlockingQueue<MSMSBlock> outputBlockQueue, SearchParameters parameters) {
		this.diaFile=diaFile;
		this.outputBlockQueue=outputBlockQueue;
		this.parameters=parameters;
		this.fraction=fraction;
	}

	@Override
	public void run() {
		try {
			library=new LibraryFile();
			library.openFile(diaFile);
			
			producePrecursors();
			produceFragments();
			
			library.close();
			
		} catch (IOException ioe) {
			throw new EncyclopediaException("DIA writing IO error!", ioe);
		} catch (SQLException sqle) {
			throw new EncyclopediaException("DIA writing SQL error!", sqle);
		} catch (DataFormatException dfe) {
			throw new EncyclopediaException("DIA writing DFE error!", dfe);
		
		}
	}

	public void producePrecursors() throws IOException, SQLException, DataFormatException {
		Connection c = library.getConnection();
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select SpectrumName, SpectrumIndex, ScanStartTime, IonInjectionTime, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray, TIC, fraction, isolationWindowLower, isolationWindowUpper from precursor");

				ArrayList<PrecursorScan> precursors=new ArrayList<PrecursorScan>();
				while (rs.next()) {
					String spectrumName=rs.getString(1);
					int spectrumIndex=rs.getInt(2);
					float scanStartTime=rs.getFloat(3);
					Float ionInjectionTime=rs.getFloat(4);
					if (rs.wasNull()) {
						ionInjectionTime=null;
					}
					int massEncodedLength=rs.getInt(5);
					double[] massArray=ByteConverter.toDoubleArray(CompressionUtils.decompress(rs.getBytes(6), massEncodedLength));
					int intensityEncodedLength=rs.getInt(7);
					float[] intensityArray=ByteConverter.toFloatArray(CompressionUtils.decompress(rs.getBytes(8), intensityEncodedLength));
					float tic=rs.getFloat(9);
					int fraction=rs.getInt(10);
					double isolationWindowLower=rs.getDouble(11);
					double isolationWindowUpper=rs.getDouble(12);

					precursors.add(new PrecursorScan(spectrumName, spectrumIndex, scanStartTime, fraction, isolationWindowLower, isolationWindowUpper, ionInjectionTime, massArray, intensityArray, tic));
					
					if (precursors.size()>=MAX_PRECURSORS_PER_BLOCK) {
						putBlock(new MSMSBlock(precursors, new ArrayList<>()));
						precursors.clear();
					}
				}

				putBlock(new MSMSBlock(precursors, new ArrayList<>()));

			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}
	
	public void produceFragments() throws IOException, SQLException, DataFormatException {
		Connection c = library.getConnection();
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select SpectrumName, PrecursorName, SpectrumIndex, ScanStartTime, IsolationWindowLower, IsolationWindowUpper, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray, IonInjectionTime, Fraction from spectra");

				final ArrayList<FragmentScan> stripes=new ArrayList<FragmentScan>();

				while (rs.next()) {
					final String spectrumName=rs.getString(1);
					final String precursorName=rs.getString(2);
					final int spectrumIndex=rs.getInt(3);
					final float scanStartTime=rs.getFloat(4);
					final float isolationWindowLower=rs.getFloat(5);
					final float isolationWindowUpper=rs.getFloat(6);
					final int massEncodedLength=rs.getInt(7);
					final byte[] massBytes=rs.getBytes(8);
					final int intensityEncodedLength=rs.getInt(9);
					final byte[] intensityBytes=rs.getBytes(10);
					Float nullableIonInjectionTime=rs.getFloat(11);
					if (rs.wasNull()) {
						nullableIonInjectionTime=null;
					}
					final Float ionInjectionTime=nullableIonInjectionTime;
					final int fraction=rs.getInt(12);
					
					double[] massArray=ByteConverter.toDoubleArray(CompressionUtils.decompress(massBytes, massEncodedLength));
					float[] intensityArray=ByteConverter.toFloatArray(CompressionUtils.decompress(intensityBytes, intensityEncodedLength));
					FragmentScan scan=new FragmentScan(spectrumName, precursorName, spectrumIndex, scanStartTime, fraction, ionInjectionTime, isolationWindowLower, isolationWindowUpper, massArray, intensityArray);
					
					stripes.add(scan);

					if (stripes.size()>=MAX_PRECURSORS_PER_BLOCK) {
						putBlock(new MSMSBlock(new ArrayList<>(), stripes));
						stripes.clear();
					}
				}
				putBlock(new MSMSBlock(new ArrayList<>(), stripes));

			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	@Override
	public void putBlock(MSMSBlock block) {
		try {
			outputBlockQueue.put(block);
		} catch (InterruptedException ie) {
			Logger.errorLine("Mzml reading interrupted!");
			Logger.errorException(ie);
		}
	}
}
