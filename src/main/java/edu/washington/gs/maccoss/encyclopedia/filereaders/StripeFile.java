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
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;

import com.google.common.util.concurrent.ThreadFactoryBuilder;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.Version;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class StripeFile extends VersionedSQLFile implements StripeFileInterface {
	private static final Version MOST_RECENT_VERSION = new Version(0, 1, 0);

	private static final String UNKNOWN_VALUE="unknown";
	public static final String FILELOCATION_ATTRIBUTE="filelocation";
	public static final String SOURCENAME_ATTRIBUTE="sourcename";
	public static final String FILENAME_ATTRIBUTE="filename";
	public static final String TOTAL_PRECURSOR_TIC_ATTRIBUTE="totalPrecursorTIC";
	public static final String GRADIENT_LENGTH_ATTRIBUTE="gradientLength";

	public static final String DIA_EXTENSION=".dia";
	
	private File userFile=null;
	private volatile String originalFileName=null;
	private File tempFile;
	private boolean isOpen=false;
	
	private final HashMap<Range, Float> ranges=new HashMap<Range, Float>();

	private final boolean isOpenFileInPlace;
	
	public StripeFile() throws IOException {
		this(false);
	}
	
	public StripeFile(boolean isOpenFileInPlace) throws IOException {
		if (!isOpenFileInPlace){
			tempFile=File.createTempFile("encyclopedia_", DIA_EXTENSION);
			tempFile.deleteOnExit();
		}
		this.isOpenFileInPlace = isOpenFileInPlace;
	}
	
	/**
	 * it's ok that this can generate races to set originalFileName, since as long as we don't overwrite with null it'll never change
	 * @return
	 * @throws IOException
	 * @throws SQLException
	 */
	public String getOriginalFileName() {
		if (originalFileName!=null) {
			return originalFileName;
		} else {
			try {
				Map<String, String> metadata=getMetadata();
				String fname=metadata.get(FILENAME_ATTRIBUTE);
				if (fname!=null) {
					Optional<String> optional=StripeFileGenerator.getBuggyFileName(fname);
					if (optional.isPresent()) {
						originalFileName=optional.get();
					} else {
						originalFileName=fname;
					}
				}
				return originalFileName;
			} catch (IOException ioe) {
				return null;
			} catch (SQLException sqle) {
				return null;
			}
		}
	}

	public CachedStripeFile cache() throws IOException, SQLException, DataFormatException {
		return cache(this);
	}

	public static CachedStripeFile cache(StripeFileInterface stripeFile) throws IOException, SQLException, DataFormatException {
		Logger.logLine("Caching precursors...");
		ArrayList<PrecursorScan> precursors=stripeFile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE);
		HashMap<Range, ArrayList<Stripe>> stripes=new HashMap<Range, ArrayList<Stripe>>();
		final Map<Range, Float> ranges = stripeFile.getRanges();
		for (Range range : ranges.keySet()) {
			Logger.logLine("Caching range "+range.toString()+"...");
			stripes.put(range, stripeFile.getStripes(range.getMiddle(), -Float.MAX_VALUE, Float.MAX_VALUE, false));
		}
		final File userFile = stripeFile.getFile();
		Logger.logLine("Finished caching "+userFile.getName());
		return new CachedStripeFile(userFile, ranges, precursors, stripes);
	}
	
	public File getFile() {
		if (userFile==null) return tempFile;
		return userFile;
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
		isOpen=true;
	}
	
	public void loadRanges() throws IOException, SQLException {
		Connection c = getConnection();
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
		Connection c = getConnection();
		try {
			PreparedStatement prep=c.prepareStatement("insert into ranges (Start, Stop, DutyCycle) VALUES (?,?,?)");
			try {
				int rangeCount=0;
				for (Entry<Range, Float> entry : ranges.entrySet()) {
					Range range=entry.getKey();
					Float value=entry.getValue();
					if (value!=null) {
						prep.setFloat(1, range.getStart());
						prep.setFloat(2, range.getStop());
						prep.setFloat(3, value);
						prep.addBatch();
						rangeCount++;
						System.out.println(value+"\t"+range);
					}
				}
				if (rangeCount>0) {
					prep.executeBatch();
				}
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
		if (isOpenFileInPlace) {
			tempFile=userFile;
		} else {
			if (userFile!=null) {
				Files.copy(userFile.toPath(), tempFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
			}
		}
		createNewTables();
	}

	public void saveAsFile(File userFile) throws IOException, SQLException {
		this.userFile=userFile;
		saveFile();
	}

	public void saveFile() throws IOException, SQLException {
		writeRanges();
		
		if (userFile!=null && !isOpenFileInPlace) {
			Files.copy(tempFile.toPath(), userFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
		}
	}

	public void setFileName(String fileName, String sourceName, String fileLocation) throws IOException, SQLException {
		HashMap<String, String> map=new HashMap<String, String>();
		map.put(FILENAME_ATTRIBUTE, fileName==null?UNKNOWN_VALUE:fileName);
		map.put(SOURCENAME_ATTRIBUTE, sourceName==null?UNKNOWN_VALUE:sourceName);
		map.put(FILELOCATION_ATTRIBUTE, fileLocation==null?UNKNOWN_VALUE:fileLocation);
		addMetadata(map);
	}
	
	public float getTIC() throws IOException, SQLException {
		String value=getMetadata().get(StripeFile.TOTAL_PRECURSOR_TIC_ATTRIBUTE);
		if (value==null) return 0.0f;
		return Float.parseFloat(value);
	}
	
	public float getGradientLength() throws IOException, SQLException {
		String value=getMetadata().get(StripeFile.GRADIENT_LENGTH_ATTRIBUTE);
		if (value==null) {
			float rt=0.0f;
			Connection c = getConnection();
			try {
				Statement s=c.createStatement();
				try {
					ResultSet rs=s.executeQuery("select max(scanstarttime) from spectra");

					while (rs.next()) {
						rt=rs.getFloat(1);
					}
				} finally {
					s.close();
				}
			} finally {
				c.close();
			}
			
			if (rt>0.0f) {
				addMetadata(StripeFile.GRADIENT_LENGTH_ATTRIBUTE, Float.toString(rt));
			}
			return rt;
		}
		return Float.parseFloat(value);
	}

	public void addPrecursor(ArrayList<PrecursorScan> precursors) throws IOException, SQLException {
		Connection c = getConnection();
		try {
			PreparedStatement prep=c.prepareStatement("insert into precursor (SpectrumName, SpectrumIndex, ScanStartTime, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray, TIC) VALUES (?,?,?,?,?,?,?, ?)");
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
					prep.setFloat(8, precursor.getTIC());
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

	protected Connection getConnection() throws IOException, SQLException {
		if (isOpenFileInPlace && !userFile.exists()){
			throw new IllegalStateException("No file to obtain a connection to!");
		}
		return isOpenFileInPlace ? getConnection(userFile): getConnection(tempFile);
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface#getPrecursors(float, float)
	 */
	@Override
	public ArrayList<PrecursorScan> getPrecursors(float minRT, float maxRT) throws IOException, SQLException,DataFormatException {
		Connection c = getConnection();
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select SpectrumName, SpectrumIndex, ScanStartTime, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray, TIC from precursor "
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
					float tic=rs.getFloat(8);

					precursors.add(new PrecursorScan(spectrumName, spectrumIndex, scanStartTime, massArray, intensityArray, tic));
				}

				return precursors;
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	public boolean isOpenFileInPlace() {
		return isOpenFileInPlace;
	}

	private static final int NUMBER_OF_STRIPES_AT_ONCE=10;
	public void addStripe(ArrayList<Stripe> stripes) throws IOException, SQLException {
		Connection c = getConnection();
		try {
			int start=0;
			int stop=NUMBER_OF_STRIPES_AT_ONCE;
			while (stop<stripes.size()) {
				internalAddStripeToConnection(stripes.subList(start, stop), c);
				start=stop;
				stop=stop+NUMBER_OF_STRIPES_AT_ONCE;
			}
			if (start<stripes.size()) {
				internalAddStripeToConnection(stripes.subList(start, stripes.size()), c);
			}

			c.commit();
			
		} finally {
			c.close();
		}
	}

	private void internalAddStripeToConnection(List<Stripe> stripes, Connection c) throws SQLException, IOException {
		StringBuilder sb=new StringBuilder("insert into spectra (SpectrumName, PrecursorName, SpectrumIndex, ScanStartTime, IsolationWindowLower, IsolationWindowCenter, IsolationWindowUpper, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray)");
		sb.append(" VALUES (?,?,?,?,?,?,?,?,?,?,?)");
		for (int i=1; i<stripes.size(); i++) {
			sb.append(", (?,?,?,?,?,?,?,?,?,?,?)");
		}
		PreparedStatement prep=c.prepareStatement(sb.toString());
		try {
			int index=1;
			for (Stripe stripe : stripes) {
				prep.setString(index++, stripe.getSpectrumName());
				prep.setString(index++, stripe.getPrecursorName());
				prep.setInt(index++, stripe.getSpectrumIndex());
				prep.setFloat(index++, stripe.getScanStartTime());
				prep.setFloat(index++, stripe.getIsolationWindowLower());
				prep.setFloat(index++, stripe.getIsolationWindowCenter());
				prep.setFloat(index++, stripe.getIsolationWindowUpper());
				byte[] massByteArray=ByteConverter.toByteArray(stripe.getMassArray());
				prep.setInt(index++, massByteArray.length);
				prep.setBytes(index++, CompressionUtils.compress(massByteArray));
				byte[] intensityByteArray=ByteConverter.toByteArray(stripe.getIntensityArray());
				prep.setInt(index++, intensityByteArray.length);
				prep.setBytes(index++, CompressionUtils.compress(intensityByteArray));
			}
			prep.execute();
			prep.close();
		} finally {
			prep.close();
		}
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface#getStripes(double, float, float, boolean)
	 */
	@Override
	public ArrayList<Stripe> getStripes(double targetMz, float minRT, float maxRT, final boolean sqrt) throws IOException, SQLException {
		Connection c = getConnection();
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select SpectrumName, PrecursorName, SpectrumIndex, ScanStartTime, IsolationWindowLower, IsolationWindowUpper, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray from spectra "
						+"where IsolationWindowLower <= "+targetMz+" and IsolationWindowUpper >= "+targetMz+" and ScanStartTime between "+minRT+" and "+maxRT);
				
				final Vector<Stripe> stripes=new Vector<Stripe>();
				
				int cores=Runtime.getRuntime().availableProcessors();
				ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("STRIPE_"+targetMz+"-%d").setDaemon(true).build();
				LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
				ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 

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
					executor.submit(new Runnable() {
						@Override
						public void run() {
							try {
								stripes.add(getStripe(sqrt, spectrumName, precursorName, spectrumIndex, scanStartTime, isolationWindowLower, isolationWindowUpper, massEncodedLength, massBytes,
										intensityEncodedLength, intensityBytes));
							} catch (DataFormatException dfe) {
								throw new EncyclopediaException(dfe);
							} catch (IOException ioe) {
								throw new EncyclopediaException(ioe);
							}
						}
					});
				}

				executor.shutdown();
				try {
					executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
				} catch (InterruptedException ie) {
					throw new EncyclopediaException(ie);
				}
				return new ArrayList<Stripe>(stripes);
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface#getStripes(Range, float, float, boolean)
	 */
	@Override
	public ArrayList<Stripe> getStripes(Range targetMzRange, float minRT, float maxRT, final boolean sqrt) throws IOException, SQLException {
		Connection c = getConnection();
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select SpectrumName, PrecursorName, SpectrumIndex, ScanStartTime, IsolationWindowLower, IsolationWindowUpper, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray from spectra "
						+"where  IsolationWindowLower <= "+targetMzRange.getStop()+" and IsolationWindowUpper >= "+targetMzRange.getStart()+" and ScanStartTime between "+minRT+" and "+maxRT);
				
				final Vector<Stripe> stripes=new Vector<Stripe>();
				
				int cores=Runtime.getRuntime().availableProcessors();
				ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("STRIPE_"+targetMzRange.getStart()+"_"+targetMzRange.getStop()+"-%d").setDaemon(true).build();
				LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
				ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 

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
					executor.submit(new Runnable() {
						@Override
						public void run() {
							try {
								stripes.add(getStripe(sqrt, spectrumName, precursorName, spectrumIndex, scanStartTime, isolationWindowLower, isolationWindowUpper, massEncodedLength, massBytes,
										intensityEncodedLength, intensityBytes));
							} catch (DataFormatException dfe) {
								throw new EncyclopediaException(dfe);
							} catch (IOException ioe) {
								throw new EncyclopediaException(ioe);
							}
						}
					});
				}

				executor.shutdown();
				try {
					executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
				} catch (InterruptedException ie) {
					throw new EncyclopediaException(ie);
				}
				return new ArrayList<Stripe>(stripes);
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	private Stripe getStripe(boolean sqrt, String spectrumName, String precursorName, int spectrumIndex, float scanStartTime, float isolationWindowLower,
			float isolationWindowUpper, int massEncodedLength, byte[] massBytes, int intensityEncodedLength, byte[] intensityBytes) throws IOException, DataFormatException {
		double[] massArray=ByteConverter.toDoubleArray(CompressionUtils.decompress(massBytes, massEncodedLength));
		float[] intensityArray=ByteConverter.toFloatArray(CompressionUtils.decompress(intensityBytes, intensityEncodedLength));
		if (sqrt) {
			intensityArray=General.protectedSqrt(intensityArray);
		}
		return new Stripe(spectrumName, precursorName, spectrumIndex, scanStartTime, isolationWindowLower, isolationWindowUpper, massArray, intensityArray);
	}

	@Override
	public Version getMostRecentVersion() {
		return MOST_RECENT_VERSION;
	}

	@Override
	protected void applyPatches(Version currentVersion, Statement s) throws IOException, SQLException {
		if (new Version(0, 0, 0).equals(currentVersion)) {
			// Because versioning was added to StripeFile later, version 0.0.0 can mean either a new file
			// or a file created before versioning. Thus we must check for the existence of tables we want to alter.
			if (doesTableExist(s.getConnection(), "precursor")) {
				s.execute("alter table precursor add column TIC float");
			}
		}
	}

	@Override
	protected void createTables(Statement s) throws IOException, SQLException {
		s.execute("create table if not exists metadata ( Key string not null, Value string not null, primary key (Key) )");
		s.execute("create table if not exists ranges ( Start float not null, Stop float not null, DutyCycle float not null )");
		s.execute("create table if not exists spectra ( SpectrumName string not null, PrecursorName string, SpectrumIndex int not null, ScanStartTime float not null, IsolationWindowLower float not null, IsolationWindowCenter float not null, IsolationWindowUpper float not null, MassEncodedLength int not null, MassArray blob not null, IntensityEncodedLength int not null, IntensityArray blob not null, primary key (SpectrumIndex) )");
		s.execute("create table if not exists precursor ( SpectrumName string not null, SpectrumIndex int not null, ScanStartTime float not null, MassEncodedLength int not null, MassArray blob not null, IntensityEncodedLength int not null, IntensityArray blob not null, TIC float, primary key (SpectrumIndex) )");

		s.execute("create index if not exists \"spectra_index_isolation_window_lower\" on \"spectra\" (\"IsolationWindowLower\" ASC)");
		s.execute("create index if not exists \"spectra_index_isolation_window_upper\" on \"spectra\" (\"IsolationWindowUpper\" ASC)");
		s.execute("create index if not exists \"spectra_index_scan_start_time\" on \"spectra\" (\"ScanStartTime\" ASC)");
		s.execute("create index if not exists \"precursor_index_scan_start_time\" on \"precursor\" (\"ScanStartTime\" ASC)");
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface#close()
	 */
	@Override
	public void close() {
		if (!isOpenFileInPlace && !tempFile.delete()) {
			Logger.errorLine("Error deleting temp file!");
		}
		isOpen=false;
	}
	@Override
	public boolean isOpen() {
		return isOpen;
	}
}
