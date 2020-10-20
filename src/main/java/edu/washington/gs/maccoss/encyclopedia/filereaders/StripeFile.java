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
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
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

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableMultimap;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.util.concurrent.ThreadFactoryBuilder;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.mzml.InstrumentComponent;
import edu.washington.gs.maccoss.encyclopedia.filereaders.mzml.InstrumentId;
import edu.washington.gs.maccoss.encyclopedia.filereaders.mzml.InstrumentMapTranscoder;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.Version;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class StripeFile extends SQLFile implements StripeFileInterface {
	public static final DateFormat m_ISO8601Local = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'Z'");
	private static final Version MOST_RECENT_VERSION = new Version(0, 4, 0);

	private static final String UNKNOWN_VALUE="unknown";
	public static final String FILELOCATION_ATTRIBUTE="filelocation";
	public static final String SOURCENAME_ATTRIBUTE="sourcename";
	public static final String FILENAME_ATTRIBUTE="filename";
	public static final String TOTAL_PRECURSOR_TIC_ATTRIBUTE="totalPrecursorTIC";
	public static final String GRADIENT_LENGTH_ATTRIBUTE="gradientLength";
	public static final String SOFTWARE_VERSION_PREFIX = "SoftwareVersion_";
	public static final String RUN_START_TIME = "runStartTime";
	public static final String SOFTWARE_VERSIONS_DELIMITER = ";";
	public static final String INSTRUMENT_CONFIGURATIONS = "InstrumentConfigurations";

	public static final String DIA_EXTENSION=".dia";

	private File userFile=null;
	private volatile String originalFileName=null;
	private File tempFile;
	private boolean isOpen=false;

	private final HashMap<Range, WindowData> ranges=new HashMap<Range, WindowData>();

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
		HashMap<Range, ArrayList<FragmentScan>> stripes=new HashMap<Range, ArrayList<FragmentScan>>();
		final Map<Range, WindowData> ranges = stripeFile.getRanges();
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
	public HashMap<Range, WindowData> getRanges() {
		return (HashMap<Range, WindowData>)ranges.clone();
	}

	public void setRanges(HashMap<Range, WindowData> ranges) {
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
				ResultSet rs=s.executeQuery("select Start, Stop, DutyCycle,NumWindows from Ranges");

				while (rs.next()) {
					float start=rs.getFloat(1);
					float stop=rs.getFloat(2);
					float dutyCycle=rs.getFloat(3);
					int numWindows=rs.getInt(4);
					ranges.put(new Range(start, stop), new WindowData(dutyCycle, numWindows));
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
			PreparedStatement prep=c.prepareStatement("insert into ranges (Start, Stop, DutyCycle, NumWindows) VALUES (?,?,?,?)");
			try {
				int rangeCount=0;
				for (Entry<Range, WindowData> entry : ranges.entrySet()) {
					Range range=entry.getKey();
					WindowData data=entry.getValue();
					if (data!=null) {
						float dutyCycle=data.getAverageDutyCycle();
						int numWindows=data.getNumberOfMSMS();
						prep.setFloat(1, range.getStart());
						prep.setFloat(2, range.getStop());
						prep.setFloat(3, dutyCycle);
						prep.setInt(4, numWindows);
						prep.addBatch();
						rangeCount++;
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

	public void setFileVersion() throws IOException, SQLException {
		HashMap<String, String> map=new HashMap<String, String>();
		map.put(VERSION_STRING, getMostRecentVersion().toString());
		addMetadata(map);
	}

	public void saveFile() throws IOException, SQLException {
		writeRanges();

		if (userFile!=null) {
			setFileVersion();

			if (!isOpenFileInPlace) {
				Files.copy(tempFile.toPath(), userFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
			}
		}
	}

	public void setFileName(String fileName, String sourceName, String fileLocation) throws IOException, SQLException {
		HashMap<String, String> map=new HashMap<String, String>();
		map.put(FILENAME_ATTRIBUTE, fileName==null?UNKNOWN_VALUE:fileName);
		map.put(SOURCENAME_ATTRIBUTE, sourceName==null?UNKNOWN_VALUE:sourceName);
		map.put(FILELOCATION_ATTRIBUTE, fileLocation==null?UNKNOWN_VALUE:fileLocation);
		addMetadata(map);
	}

	public void setInstrumentConfiguration(ImmutableMultimap<InstrumentId, InstrumentComponent> instrumentConfigurations) throws IOException, SQLException {
		if (!instrumentConfigurations.isEmpty()) {
			Map<String, String> m = Maps.newHashMap();
			m.put(INSTRUMENT_CONFIGURATIONS, InstrumentMapTranscoder.encode(instrumentConfigurations));
			addMetadata(m);
		}
	}
	
	public void setStartTime(Date startTime) throws IOException, SQLException {
		addMetadata(RUN_START_TIME, m_ISO8601Local.format(startTime));
	}

	public void setSoftwareVersions(final Multimap<String, String> softwareAccessionIdToVersion) throws IOException, SQLException {
		if (!softwareAccessionIdToVersion.isEmpty()) {
			Map<String, String> toAdd = Maps.newHashMap();
			softwareAccessionIdToVersion.asMap().forEach((key, value) -> {
				toAdd.put(SOFTWARE_VERSION_PREFIX + key, Joiner.on(SOFTWARE_VERSIONS_DELIMITER).join(value));
			});
			addMetadata(toAdd);
		}
	}

	public void addMetadata(String key, String value) throws IOException, SQLException {
		HashMap<String, String> map=new HashMap<String, String>();
		map.put(key, value==null?UNKNOWN_VALUE:value);
		addMetadata(map);
	}

	public HashMap<String, String> getMetadata() throws IOException, SQLException {
		Connection c = getConnection();
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

	public void addMetadata(Map<String, String> data) throws IOException, SQLException {
		Connection c = getConnection();
		try {
			PreparedStatement prep=c.prepareStatement("insert or replace into metadata (Key, Value) VALUES (?,?)");
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
		Connection c = getConnection();
		try {
			PreparedStatement prep=c.prepareStatement("insert into precursor (SpectrumName, SpectrumIndex, ScanStartTime, IonInjectionTime, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray, TIC, Fraction, IsolationWindowLower, IsolationWindowUpper) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)");
			try {
				for (PrecursorScan precursor : precursors) {
					prep.setString(1, precursor.getSpectrumName());
					prep.setInt(2, precursor.getSpectrumIndex());
					prep.setFloat(3, precursor.getScanStartTime());
					prep.setFloat(4, precursor.getIonInjectionTime());
					byte[] massByteArray=ByteConverter.toByteArray(precursor.getMassArray());
					prep.setInt(5, massByteArray.length);
					prep.setBytes(6, CompressionUtils.compress(massByteArray));
					byte[] intensityByteArray=ByteConverter.toByteArray(precursor.getIntensityArray());
					prep.setInt(7, intensityByteArray.length);
					prep.setBytes(8, CompressionUtils.compress(intensityByteArray));
					prep.setFloat(9, precursor.getTIC());
					prep.setInt(10, precursor.getFraction());
					prep.setDouble(11, precursor.getIsolationWindowLower());
					prep.setDouble(12, precursor.getIsolationWindowUpper());
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

	public Connection getConnection() throws IOException, SQLException {
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
				ResultSet rs=s.executeQuery("select SpectrumName, SpectrumIndex, ScanStartTime, IonInjectionTime, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray, TIC, fraction, isolationWindowLower, isolationWindowUpper from precursor "
						+"where ScanStartTime between "+minRT+" and "+maxRT);

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
	public void addStripe(ArrayList<FragmentScan> stripes) throws IOException, SQLException {
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

	private void internalAddStripeToConnection(List<FragmentScan> stripes, Connection c) throws SQLException, IOException {
		StringBuilder sb=new StringBuilder("insert into spectra (SpectrumName, PrecursorName, SpectrumIndex, ScanStartTime, Fraction, IonInjectionTime, IsolationWindowLower, IsolationWindowCenter, IsolationWindowUpper, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray)");
		sb.append(" VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)");
		for (int i=1; i<stripes.size(); i++) {
			sb.append(", (?,?,?,?,?,?,?,?,?,?,?,?,?)");
		}
		PreparedStatement prep=c.prepareStatement(sb.toString());
		try {
			int index=1;
			for (FragmentScan stripe : stripes) {
				prep.setString(index++, stripe.getSpectrumName());
				prep.setString(index++, stripe.getPrecursorName());
				prep.setInt(index++, stripe.getSpectrumIndex());
				prep.setFloat(index++, stripe.getScanStartTime());
				prep.setInt(index++, stripe.getFraction());
				prep.setFloat(index++, stripe.getIonInjectionTime());
				prep.setDouble(index++, stripe.getIsolationWindowLower());
				prep.setDouble(index++, stripe.getIsolationWindowCenter());
				prep.setDouble(index++, stripe.getIsolationWindowUpper());
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
	public ArrayList<FragmentScan> getStripes(double targetMz, float minRT, float maxRT, final boolean sqrt) throws IOException, SQLException {
		Connection c = getConnection();
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select SpectrumName, PrecursorName, SpectrumIndex, ScanStartTime, IsolationWindowLower, IsolationWindowUpper, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray, IonInjectionTime, Fraction from spectra "
						+"where IsolationWindowLower <= "+targetMz+" and IsolationWindowUpper >= "+targetMz+" and ScanStartTime between "+minRT+" and "+maxRT);

				final Vector<FragmentScan> stripes=new Vector<FragmentScan>();

				int cores=Runtime.getRuntime().availableProcessors();
				ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("STRIPE_"+targetMz+"-%d").setDaemon(true).build();
				LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
				ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory);

				while (rs.next()) {
					final String spectrumName=rs.getString(1);
					final String precursorName=rs.getString(2);
					final int spectrumIndex=rs.getInt(3);
					final float scanStartTime=rs.getFloat(4);
					final double isolationWindowLower=rs.getDouble(5);
					final double isolationWindowUpper=rs.getDouble(6);
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
					executor.submit(new Runnable() {
						@Override
						public void run() {
							try {
								stripes.add(getStripe(sqrt, spectrumName, precursorName, spectrumIndex, scanStartTime, fraction, ionInjectionTime, isolationWindowLower, isolationWindowUpper, massEncodedLength, massBytes,
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
				return new ArrayList<FragmentScan>(stripes);
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
	public ArrayList<FragmentScan> getStripes(Range targetMzRange, float minRT, float maxRT, final boolean sqrt) throws IOException, SQLException {
		Connection c = getConnection();
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select SpectrumName, PrecursorName, SpectrumIndex, ScanStartTime, IsolationWindowLower, IsolationWindowUpper, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray, IonInjectionTime, Fraction from spectra "
						+"where  IsolationWindowLower <= "+targetMzRange.getStop()+" and IsolationWindowUpper >= "+targetMzRange.getStart()+" and ScanStartTime between "+minRT+" and "+maxRT);

				final Vector<FragmentScan> stripes=new Vector<FragmentScan>();

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
					Float nullableIonInjectionTime=rs.getFloat(11);
					if (rs.wasNull()) {
						nullableIonInjectionTime=null;
					}
					final Float ionInjectionTime=nullableIonInjectionTime;
					final int fraction=rs.getInt(12);
					
					executor.submit(new Runnable() {
						@Override
						public void run() {
							try {
								stripes.add(getStripe(sqrt, spectrumName, precursorName, spectrumIndex, scanStartTime, fraction, ionInjectionTime, isolationWindowLower, isolationWindowUpper, massEncodedLength, massBytes,
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
				return new ArrayList<FragmentScan>(stripes);
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	private FragmentScan getStripe(boolean sqrt, String spectrumName, String precursorName, int spectrumIndex, Float scanStartTime, int fraction, Float ionInjectionTime, double isolationWindowLower,
			double isolationWindowUpper, int massEncodedLength, byte[] massBytes, int intensityEncodedLength, byte[] intensityBytes) throws IOException, DataFormatException {
		double[] massArray=ByteConverter.toDoubleArray(CompressionUtils.decompress(massBytes, massEncodedLength));
		float[] intensityArray=ByteConverter.toFloatArray(CompressionUtils.decompress(intensityBytes, intensityEncodedLength));
		if (sqrt) {
			intensityArray=General.protectedSqrt(intensityArray);
		}
		return new FragmentScan(spectrumName, precursorName, spectrumIndex, scanStartTime, fraction, ionInjectionTime, isolationWindowLower, isolationWindowUpper, massArray, intensityArray);
	}

	public Version getVersion() throws IOException, SQLException {
		HashMap<String, String> meta=getMetadata();
		return new Version(meta.get(VERSION_STRING));
	}

	public Version getMostRecentVersion() {
		return MOST_RECENT_VERSION;
	}

	protected void applyPatches(Version currentVersion, Statement s) throws IOException, SQLException {

			if (new Version(0, 1, 0).amIAbove(currentVersion)) {
				s.execute("alter table precursor add column TIC float");
				s.getConnection().commit();
				populateTICColumn(s.getConnection());
			}
			if (new Version(0, 2, 0).amIAbove(currentVersion)) {
				s.execute("alter table precursor add column IonInjectionTime float");
				s.execute("alter table spectra add column IonInjectionTime float");
				s.getConnection().commit();
			}
			if (new Version(0, 3, 0).amIAbove(currentVersion)) {
				s.execute("alter table spectra add column fraction int");
				s.execute("alter table precursor add column fraction int");
				s.execute("alter table precursor add column IsolationWindowLower float");
				s.execute("alter table precursor add column IsolationWindowUpper float");
				s.execute("update spectra set fraction=0");
				s.execute("update precursor set fraction=0,IsolationWindowLower=0,IsolationWindowUpper=999999999");
				s.getConnection().commit();
			}
			if (new Version(0, 4, 0).amIAbove(currentVersion)) {
				s.execute("alter table ranges add column numWindows int");
			}

	}

	private void populateTICColumn(Connection connection) throws SQLException, IOException {
		final boolean wasAutoCommit = connection.getAutoCommit();
		connection.setAutoCommit(false);
		Statement s1=null;
		PreparedStatement batchInsertStatement=null;
		try {
			s1 = connection.createStatement();
			s1.execute("create table tic_temp_store (SpectrumIndex int primary key, TIC float)");
			batchInsertStatement = connection.prepareStatement("insert into tic_temp_store values (?, ?)");
			final ResultSet resultSet = s1.executeQuery("select SpectrumIndex, IntensityEncodedLength, IntensityArray from precursor");

			int count=0;

			while (resultSet.next()) {
				try {
					final int spectrumIndex = resultSet.getInt(1);
					final int intensityEncodedLength = resultSet.getInt(2);
					final float[] intensityArray = ByteConverter.toFloatArray(CompressionUtils.decompress(
														resultSet.getBytes(3), intensityEncodedLength));
					final float tic = General.sum(intensityArray);

					batchInsertStatement.setInt(1, spectrumIndex);
					batchInsertStatement.setFloat(2, tic);
					batchInsertStatement.addBatch();

					if (++count % 8192 == 0) {
						batchInsertStatement.executeBatch();
					}
				} catch (DataFormatException e) {
					Logger.errorException(e); // log and continue
				}
			}
			batchInsertStatement.executeBatch();

			connection.commit();
			s1.execute("update precursor set TIC = (select TIC from tic_temp_store where precursor.SpectrumIndex = tic_temp_store.SpectrumIndex)");

			s1.execute("drop table tic_temp_store");
			connection.commit();
		} finally {
			connection.setAutoCommit(wasAutoCommit);
			if (s1!=null) {
				s1.close();
			}
			if (batchInsertStatement!=null) {
				batchInsertStatement.close();
			}
		}
	}

	private void createNewTables() throws IOException, SQLException {
		Connection c = getConnection();
		try {
			Statement s=c.createStatement();
			try {
				Version version = doesTableExist(c, "metadata") ? getVersion() : null;
				if (version!=null) {
					applyPatches(version, s);
				}

				s.execute("create table if not exists metadata ( Key string not null, Value string not null, primary key (Key) )");
				s.execute("create table if not exists ranges ( Start float not null, Stop float not null, DutyCycle float not null, NumWindows int )");
				s.execute("create table if not exists spectra ( Fraction int not null, SpectrumName string not null, PrecursorName string, SpectrumIndex int not null, ScanStartTime float not null, IonInjectionTime float, IsolationWindowLower float not null, IsolationWindowCenter float not null, IsolationWindowUpper float not null, MassEncodedLength int not null, MassArray blob not null, IntensityEncodedLength int not null, IntensityArray blob not null, primary key (SpectrumIndex) )");
				s.execute("create table if not exists precursor ( Fraction int not null, SpectrumName string not null, SpectrumIndex int not null, ScanStartTime float not null, IonInjectionTime float, IsolationWindowLower float not null, IsolationWindowUpper float not null, MassEncodedLength int not null, MassArray blob not null, IntensityEncodedLength int not null, IntensityArray blob not null, TIC float, primary key (SpectrumIndex) )");

				s.execute("create index if not exists \"spectra_index_isolation_window_lower\" on \"spectra\" (\"IsolationWindowLower\" ASC)");
				s.execute("create index if not exists \"spectra_index_isolation_window_upper\" on \"spectra\" (\"IsolationWindowUpper\" ASC)");
				s.execute("create index if not exists \"spectra_index_scan_start_time\" on \"spectra\" (\"ScanStartTime\" ASC)");

				s.execute("create index if not exists \"precursor_index_isolation_window_lower\" on \"precursor\" (\"IsolationWindowLower\" ASC)");
				s.execute("create index if not exists \"precursor_index_isolation_window_upper\" on \"precursor\" (\"IsolationWindowUpper\" ASC)");
				s.execute("create index if not exists \"precursor_index_scan_start_time\" on \"precursor\" (\"ScanStartTime\" ASC)");

				c.commit();
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
		if (isOpenFileInPlace) {
			setFileVersion();
		}
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