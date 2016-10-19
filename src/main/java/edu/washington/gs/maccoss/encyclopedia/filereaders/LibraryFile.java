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
import java.sql.Types;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.StringTokenizer;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeakLocationInferrer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Chromatogram;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ChromatogramLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.io.Version;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class LibraryFile extends SQLFile implements LibraryInterface {
	private static final String SOURCEFILE_TIC_PREFIX="TIC_";
	private static final String SOURCEFILE_STRING="sourcefile";
	private static final String SOURCE_FILE_SPLIT="|";
	public static final String ELIB=".elib";
	public static final String VERSION_STRING="version";
	public static final Version[] ACCEPTABLE_VERSIONS=new Version[] {new Version(0, 1, 0), new Version(0, 1, 1), new Version(0, 1, 2), new Version(0, 1, 3), new Version(0, 1, 4), new Version(0, 1, 5)};
	public static final Version MOST_RECENT_VERSION=new Version(0, 1, 5);

	private File userFile=null;
	private final File tempFile;

	public LibraryFile() throws IOException {
		tempFile=File.createTempFile("encyclopedia_", ELIB);
		tempFile.deleteOnExit();
	}

	public static boolean isVersionAcceptable(Version version) {
		if (version==null) return false;

		for (Version string : ACCEPTABLE_VERSIONS) {
			if (string.equals(version)) return true;
		}
		return false;
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
	
	public File getFile() {
		return userFile;
	}

	public void saveFile() throws IOException, SQLException {
		if (userFile!=null) {
			setFileVersion();

			Connection c=getConnection();

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

	public void addTIC(StripeFileInterface diaFile) throws IOException, SQLException {
		String key=SOURCEFILE_TIC_PREFIX+diaFile.getOriginalFileName();

		HashMap<String, String> map=new HashMap<String, String>();
		map.put(key, Float.toString(diaFile.getTIC()));

		addMetadata(map);
	}

	public float getTIC(StripeFileInterface diaFile) throws IOException, SQLException {
		String originalFileName=diaFile.getOriginalFileName();
		return getTIC(originalFileName);
	}

	public float getTIC(String originalFileName) throws IOException, SQLException {
		String key=SOURCEFILE_TIC_PREFIX+originalFileName;

		String value=getMetadata().get(key);
		if (value==null) return 0.0f;
		return Float.parseFloat(value);
	}

	public void setSources(ArrayList<SearchJobData> sources) throws IOException, SQLException {
		HashMap<String, String> map=new HashMap<String, String>();
		StringBuilder sb=new StringBuilder();
		for (SearchJobData searchJobData : sources) {
			if (sb.length()>0) {
				sb.append(SOURCE_FILE_SPLIT);
			}
			sb.append(searchJobData.getDiaFile().getAbsolutePath());

		}
		map.put(SOURCEFILE_STRING, sb.toString());
		addMetadata(map);
	}

	public void setFileVersion() throws IOException, SQLException {
		HashMap<String, String> map=new HashMap<String, String>();
		map.put(VERSION_STRING, MOST_RECENT_VERSION.toString());
		addMetadata(map);
	}

	public void addMetadata(Map<String, String> data) throws IOException, SQLException {
		Connection c=getConnection();
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

	public Version getVersion() throws IOException, SQLException {
		HashMap<String, String> meta=getMetadata();
		return new Version(meta.get(VERSION_STRING));
	}

	public ArrayList<File> getSourceFiles() throws IOException, SQLException {
		HashMap<String, String> meta=getMetadata();
		String sources=meta.get(SOURCEFILE_STRING);
		if (sources==null) return new ArrayList<File>();

		StringTokenizer st=new StringTokenizer(sources, SOURCE_FILE_SPLIT);
		ArrayList<File> files=new ArrayList<File>();
		while (st.hasMoreTokens()) {
			files.add(new File(st.nextToken()));
		}
		return files;
	}

	public Optional<StripeFileInterface> getSource(SearchParameters parameters) {
		try {
			ArrayList<File> files=getSourceFiles();
			if (files.size()==0||files.size()>1) return Optional.empty();

			StripeFileInterface file=MzmlToDIAConverter.getFile(files.get(0), parameters);
			return Optional.ofNullable(file);
		} catch (Exception e) {
			return Optional.empty();
		}
	}

	public HashMap<String, String> getMetadata() throws IOException, SQLException {
		Connection c=getConnection();
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

	public void addIntegratedEntries(ArrayList<IntegratedLibraryEntry> entries, Optional<PeakLocationInferrer> inferrer) throws IOException, SQLException {
		// first add normal data
		HashMap<String, LibraryEntry> repeatsCatcher=new HashMap<String, LibraryEntry>();
		for (IntegratedLibraryEntry entry : entries) {
			String key=entry.getPeptideModSeq()+"+"+entry.getPrecursorCharge()+","+entry.getSource();
			LibraryEntry prev=repeatsCatcher.get(key);
			if (prev==null) {
				repeatsCatcher.put(key, entry);
			} else {
				Logger.errorLine("Found collision writing elib: "+key+" ("+entry.getScore()+" vs"+prev.getScore()+"), keeping best scoring. Let Brian know if you see this!");
				if (entry.getScore()>prev.getScore()) {
					repeatsCatcher.put(key, entry);
				}
			}
		}
		ArrayList<LibraryEntry> uniqueEntries=new ArrayList<LibraryEntry>(repeatsCatcher.values());
		addEntries(uniqueEntries);

		// then add integrated areas
		Connection c=getConnection();
		try {
			PreparedStatement peptidePrep=c.prepareStatement(
					"INSERT INTO peptidequants (PrecursorCharge, PeptideModSeq, PeptideSeq, SourceFile, LocalizationPeptideModSeq, LocalizationScore, NumberOfMods, IsSiteSpecific, RTInSecondsCenter, RTInSecondsStart, RTInSecondsStop, TotalIntensity, NumberOfQuantIons, BestFragmentCorrelation, BestFragmentDeltaMassPPM, MedianChromatogramEncodedLength, MedianChromatogramArray) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
			PreparedStatement fragmentPrep=c.prepareStatement(
					"INSERT INTO fragmentquants (PrecursorCharge, PeptideModSeq, PeptideSeq, SourceFile, IonType, FragmentMass, Correlation, Background, DeltaMassPPM, Intensity) VALUES (?,?,?,?,?,?,?,?,?,?)");
			try {
				for (LibraryEntry recast : uniqueEntries) {
					IntegratedLibraryEntry entry=(IntegratedLibraryEntry)recast;
					String sourceFile=entry.getSource();
					TransitionRefinementData data=entry.getRefinementData();
					
					HashMap<String, TransitionRefinementData> uniqueDataMap=new HashMap<String, TransitionRefinementData>();

					if (data.getModificationQuantData().isPresent()) {
						HashMap<String, TransitionRefinementData> forms=data.getModificationQuantData().get();
						// preserve only the least ambiguous forms
						for (Entry<String, TransitionRefinementData> mapEntry : forms.entrySet()) {
							TransitionRefinementData value=mapEntry.getValue();
							
							if (uniqueDataMap.containsKey(value.getPeptideModSeq())) {
								Optional<ModificationLocalizationData> prevLocalizationData=uniqueDataMap.get(value.getPeptideModSeq()).getLocalizationData();
								int prevAmbiguityScore=prevLocalizationData.isPresent()?prevLocalizationData.get().getLocalizationPeptideModSeq().numAmbigousResidues():0;
								int newAmbiguityScore=value.getLocalizationData().isPresent()?value.getLocalizationData().get().getLocalizationPeptideModSeq().numAmbigousResidues():0;
								if (newAmbiguityScore<prevAmbiguityScore) {
									// new is less ambiguous
									uniqueDataMap.put(value.getPeptideModSeq(), value);
								}
							} else {
								uniqueDataMap.put(value.getPeptideModSeq(), value);
							}
						}
						if (forms.size()==0) {
							// always replace with perfect forms
							uniqueDataMap.put(data.getPeptideModSeq(), data);
						}
					} else {
						// always replace with perfect forms
						uniqueDataMap.put(data.getPeptideModSeq(), data);
					}
					
					for (TransitionRefinementData uniqueData : uniqueDataMap.values()) {
						prepareQuantData(uniqueData, sourceFile, inferrer, peptidePrep, fragmentPrep);
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

	private static final HashMap<String, String> equality=new HashMap<String, String>();
	public void prepareQuantData(TransitionRefinementData data, String sourceFile, Optional<PeakLocationInferrer> inferrer, PreparedStatement peptidePrep, PreparedStatement fragmentPrep)
			throws SQLException, IOException {
		float[] correlationArray=data.getCorrelationArray();
		float[] integrationArray=data.getIntegrationArray();
		float[] backgroundArray=data.getBackgroundArray();

		double[] fragmentMassArray=data.getFragmentMassArray();
		float[] deltaMassArray=data.getDeltaMassArray().get();
		float[] ppmArray=new float[deltaMassArray.length];

		float bestCorrelation=-1.0f;
		float bestDeltaMass=10.0f;
		for (int i=0; i<deltaMassArray.length; i++) {
			ppmArray[i]=deltaMassArray[i]*1000000.0f/(float)fragmentMassArray[i];
			if (correlationArray[i]>bestCorrelation) {
				bestCorrelation=correlationArray[i];
				bestDeltaMass=ppmArray[i];
			}
		}

		Pair<Float, Integer> topN;
		if (inferrer.isPresent()) {
			topN=inferrer.get().getTopNIntensity(data);
		} else {
			topN=data.getTopNIntensity(TransitionRefiner.quantitativeCorrelationThreshold, Integer.MAX_VALUE);
		}
		
		if (equality.containsKey(data.getPeptideModSeq())) {
			System.out.println("FOUND EXTERNAL COLLISION! "+data.getPeptideModSeq());
			System.out.println("PREV: "+equality.get(data.getPeptideModSeq()));
			if (data.getLocalizationData().isPresent()) {
				ModificationLocalizationData modData=data.getLocalizationData().get();
				System.out.println("NEW:  "+modData.getLocalizationPeptideModSeq().getPeptideAnnotation());
			} else {
				System.out.println("NEW:  NO LOC: "+data.getPeptideModSeq());
			}
			System.exit(1);

		} else {
			if (data.getLocalizationData().isPresent()) {
				ModificationLocalizationData modData=data.getLocalizationData().get();
				equality.put(data.getPeptideModSeq(), modData.getLocalizationPeptideModSeq().getPeptideAnnotation());
			} else {
				equality.put(data.getPeptideModSeq(), "NO LOC: "+data.getPeptideModSeq());
			}
		}

		peptidePrep.setInt(1, data.getPrecursorCharge());
		peptidePrep.setString(2, data.getPeptideModSeq());
		peptidePrep.setString(3, data.getPeptideSeq());
		peptidePrep.setString(4, sourceFile);
		
		if (data.getLocalizationData().isPresent()) {
			ModificationLocalizationData modData=data.getLocalizationData().get();
			peptidePrep.setString(5, modData.getLocalizationPeptideModSeq().getPeptideAnnotation());
			peptidePrep.setFloat(6, modData.getLocalizationScore());
			peptidePrep.setInt(7, modData.getNumberOfMods());
			peptidePrep.setBoolean(8, modData.isSiteSpecific());
			peptidePrep.setFloat(9,  modData.getRetentionTimeApexInSeconds());
			
		} else {
			peptidePrep.setNull(5, Types.VARCHAR);
			peptidePrep.setNull(6, Types.FLOAT);
			peptidePrep.setNull(7, Types.INTEGER);
			peptidePrep.setNull(8, Types.BOOLEAN);
			peptidePrep.setFloat(9,  data.getApexRT());
		}
		
		peptidePrep.setFloat(10, data.getRange().getStart());
		peptidePrep.setFloat(11, data.getRange().getStop());
		peptidePrep.setFloat(12, topN.x);
		peptidePrep.setInt(13, topN.y);
		peptidePrep.setFloat(14, bestCorrelation);
		peptidePrep.setFloat(15, bestDeltaMass);
		byte[] intensityByteArray=ByteConverter.toByteArray(data.getMedianChromatogram());
		peptidePrep.setInt(16, intensityByteArray.length);
		peptidePrep.setBytes(17, CompressionUtils.compress(intensityByteArray));
		peptidePrep.addBatch();

		for (int i=0; i<correlationArray.length; i++) {
			if (correlationArray[i]>=TransitionRefiner.identificationCorrelationThreshold) {
				fragmentPrep.setInt(1, data.getPrecursorCharge());
				fragmentPrep.setString(2, data.getPeptideModSeq());
				fragmentPrep.setString(3, data.getPeptideSeq());
				fragmentPrep.setString(4, sourceFile);
				fragmentPrep.setString(5, "BY");
				fragmentPrep.setDouble(6, fragmentMassArray[i]);
				fragmentPrep.setFloat(7, correlationArray[i]);
				fragmentPrep.setFloat(8, backgroundArray[i]);
				fragmentPrep.setFloat(9, ppmArray[i]);
				fragmentPrep.setFloat(10, integrationArray[i]);
				fragmentPrep.addBatch();
			}
		}
	}

	public void addEntries(ArrayList<LibraryEntry> entries) throws IOException, SQLException {
		Connection c=getConnection();
		try {
			PreparedStatement prep=c.prepareStatement(
					"INSERT INTO entries (PrecursorMZ, PrecursorCharge, PeptideModSeq, PeptideSeq, Copies, RTInSeconds, Score, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray, CorrelationEncodedLength, CorrelationArray, RTInSecondsStart, RTInSecondsStop, MedianChromatogramEncodedLength, MedianChromatogramArray, SourceFile) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
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

					if (entry.getMassArray().length!=entry.getIntensityArray().length)
						Logger.errorLine("MASS/INTENSITY EQUATION WRITE ERROR! "+entry.getMassArray().length+" != "+entry.getIntensityArray().length+" FOR "+entry.getPeptideModSeq()); // FIXME
					assert (entry.getMassArray().length==entry.getIntensityArray().length);

					if (entry instanceof Chromatogram) {
						Chromatogram cast=(Chromatogram)entry;

						byte[] correlationByteArray=ByteConverter.toByteArray(cast.getCorrelationArray());
						if (entry.getMassArray().length!=cast.getCorrelationArray().length)
							Logger.errorLine("MASS/CORRELATION EQUATION WRITE ERROR! "+entry.getMassArray().length+" != "+cast.getCorrelationArray().length+" FOR "+entry.getPeptideModSeq()); // FIXME
						assert (entry.getMassArray().length==cast.getCorrelationArray().length);

						prep.setInt(12, correlationByteArray.length);
						prep.setBytes(13, CompressionUtils.compress(correlationByteArray));
						prep.setFloat(14, cast.getRtRange().getStart());
						prep.setFloat(15, cast.getRtRange().getStop());

						byte[] chromatogramByteArray=ByteConverter.toByteArray(cast.getMedianChromatogram());
						prep.setInt(16, chromatogramByteArray.length);
						prep.setBytes(17, CompressionUtils.compress(chromatogramByteArray));
					} else {
						prep.setNull(12, Types.INTEGER);
						prep.setNull(13, Types.BLOB);
						prep.setNull(14, Types.FLOAT);
						prep.setNull(15, Types.FLOAT);
						prep.setNull(16, Types.INTEGER);
						prep.setNull(17, Types.BLOB);
					}

					prep.setString(18, entry.getSource());
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

	public Connection getConnection() throws IOException {
		return getConnection(tempFile);
	}

	public HashMap<PeptidePrecursor, ArrayList<LibraryEntry>> getEntries(ArrayList<PeptidePrecursor> entries, boolean sqrt) throws IOException, SQLException, DataFormatException {
		HashMap<PeptidePrecursor, ArrayList<LibraryEntry>> map=new HashMap<PeptidePrecursor, ArrayList<LibraryEntry>>();

		Connection c=getConnection();
		try {
			PreparedStatement prep=c.prepareStatement(
					"select e.PrecursorMZ, e.PrecursorCharge, e.PeptideModSeq, e.Copies, e.RTInSeconds, e.Score, e.MassEncodedLength, e.MassArray, e.IntensityEncodedLength, e.IntensityArray, e.CorrelationEncodedLength, e.CorrelationArray blob, e.RTInSecondsStart, e.RTInSecondsStop, e.MedianChromatogramEncodedLength, e.MedianChromatogramArray, p.ProteinAccessions, e.SourceFile from entries e, proteins p"
							+" where e.PeptideSeq=p.PeptideSeq and e.PeptideModSeq = ? and e.PrecursorCharge = ?");
			try {
				for (PeptidePrecursor precursor : entries) {
					prep.setString(1, precursor.getPeptideModSeq());
					prep.setByte(2, precursor.getPrecursorCharge());

					ResultSet rs=prep.executeQuery();
					ArrayList<LibraryEntry> entry=extractEntries(sqrt, rs);
					map.put(precursor, entry);
				}
			} finally {
				prep.close();
			}
		} finally {
			c.close();
		}

		return map;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface#
	 * getEntries(edu.washington.gs.maccoss.encyclopedia.datastructures.Range)
	 */
	@Override
	public ArrayList<LibraryEntry> getEntries(String peptideModSeq, byte charge, boolean sqrt) throws IOException, SQLException, DataFormatException {
		Connection c=getConnection();
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery(
						"select e.PrecursorMZ, e.PrecursorCharge, e.PeptideModSeq, e.Copies, e.RTInSeconds, e.Score, e.MassEncodedLength, e.MassArray, e.IntensityEncodedLength, e.IntensityArray, e.CorrelationEncodedLength, e.CorrelationArray blob, e.RTInSecondsStart, e.RTInSecondsStop, e.MedianChromatogramEncodedLength, e.MedianChromatogramArray, p.ProteinAccessions, e.SourceFile from entries e, proteins p"
								+" where e.PeptideSeq=p.PeptideSeq and e.PeptideModSeq = \""+peptideModSeq+"\" and e.PrecursorCharge = "+charge);

				ArrayList<LibraryEntry> entry=extractEntries(sqrt, rs);

				return entry;
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	private ArrayList<LibraryEntry> extractEntries(boolean sqrt, ResultSet rs) throws SQLException, IOException, DataFormatException {
		String peptideModSeq;
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

			if (massArray.length!=intensityArray.length) Logger.errorLine("MASS/INTENSITY EQUATION READ ERROR "+massArray.length+" != "+intensityArray.length+" FOR "+peptideModSeq); // FIXME
			assert (massArray.length==intensityArray.length);

			float[] correlationArray;
			float rtInSecondsStart;
			float rtInSecondsStop;
			float[] medianChromatogramArray;

			int correlationEncodedLength=rs.getInt(11);
			if (correlationEncodedLength==0) {
				// 0 indicates null, which indicates missing
				correlationArray=null;
				rtInSecondsStart=0.0f;
				rtInSecondsStop=0.0f;
				medianChromatogramArray=null;
			} else {
				correlationArray=ByteConverter.toFloatArray(CompressionUtils.decompress(rs.getBytes(12), correlationEncodedLength));
				if (massArray.length!=correlationArray.length) Logger.errorLine("MASS/CORRELATION EQUATION READ ERROR! "+massArray.length+" != "+correlationArray.length+" FOR "+peptideModSeq); // FIXME
				assert (massArray.length==correlationArray.length);

				rtInSecondsStart=rs.getFloat(13);
				rtInSecondsStop=rs.getFloat(14);
				int medianChromatogramEncodedLength=rs.getInt(15);
				medianChromatogramArray=ByteConverter.toFloatArray(CompressionUtils.decompress(rs.getBytes(16), medianChromatogramEncodedLength));
			}

			HashSet<String> accessions=PSMData.stringToAccessions(rs.getString(17));
			String sourceFile=rs.getString(18);
			if (correlationEncodedLength==0) {
				entry.add(new LibraryEntry(sourceFile, accessions, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray));
			} else {
				entry.add(new ChromatogramLibraryEntry(sourceFile, accessions, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray,
						correlationArray, medianChromatogramArray, new Range(rtInSecondsStart, rtInSecondsStop)));
			}
		}
		return entry;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface#
	 * getEntries(edu.washington.gs.maccoss.encyclopedia.datastructures.Range)
	 */
	@Override
	public Range getMinMaxMZ() throws IOException, SQLException {
		Connection c=getConnection();
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface#
	 * getEntries(edu.washington.gs.maccoss.encyclopedia.datastructures.Range)
	 */
	@Override
	public ArrayList<LibraryEntry> getEntries(Range precursorMz, boolean sqrt) throws IOException, SQLException, DataFormatException {
		Connection c=getConnection();
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery(
						"select e.PrecursorMZ, e.PrecursorCharge, e.PeptideModSeq, e.Copies, e.RTInSeconds, e.Score, e.MassEncodedLength, e.MassArray, e.IntensityEncodedLength, e.IntensityArray, e.CorrelationEncodedLength, e.CorrelationArray blob, e.RTInSecondsStart, e.RTInSecondsStop, e.MedianChromatogramEncodedLength, e.MedianChromatogramArray, p.ProteinAccessions, e.SourceFile from entries e, proteins p"
								+" where e.PeptideSeq=p.PeptideSeq and e.PrecursorMz between "+precursorMz.getStart()+" and "+precursorMz.getStop());

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

					float[] correlationArray;
					float rtInSecondsStart;
					float rtInSecondsStop;
					float[] medianChromatogramArray;

					int correlationEncodedLength=rs.getInt(11);
					if (correlationEncodedLength==0) {
						// 0 indicates null, which indicates missing
						correlationArray=null;
						rtInSecondsStart=0.0f;
						rtInSecondsStop=0.0f;
						medianChromatogramArray=null;
					} else {
						correlationArray=ByteConverter.toFloatArray(CompressionUtils.decompress(rs.getBytes(12), correlationEncodedLength));
						rtInSecondsStart=rs.getFloat(13);
						rtInSecondsStop=rs.getFloat(14);
						int medianChromatogramEncodedLength=rs.getInt(15);
						medianChromatogramArray=ByteConverter.toFloatArray(CompressionUtils.decompress(rs.getBytes(16), medianChromatogramEncodedLength));
					}

					String proteinToken=rs.getString(17);
					HashSet<String> accessions=PSMData.stringToAccessions(proteinToken);
					String sourceFile=rs.getString(18);
					if (correlationEncodedLength==0) {
						entry.add(new LibraryEntry(sourceFile, accessions, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray));
					} else {
						entry.add(new ChromatogramLibraryEntry(sourceFile, accessions, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray,
								correlationArray, medianChromatogramArray, new Range(rtInSecondsStart, rtInSecondsStop)));
					}
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
		Connection c=getConnection();
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery(
						"select e.PrecursorMZ, e.PrecursorCharge, e.PeptideModSeq, e.Copies, e.RTInSeconds, e.Score, e.MassEncodedLength, e.MassArray, e.IntensityEncodedLength, e.IntensityArray, e.CorrelationEncodedLength, e.CorrelationArray blob, e.RTInSecondsStart, e.RTInSecondsStop, e.MedianChromatogramEncodedLength, e.MedianChromatogramArray, p.ProteinAccessions, e.SourceFile from entries e, proteins p"
								+" where e.PeptideSeq=p.PeptideSeq");

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

					float[] correlationArray;
					float rtInSecondsStart;
					float rtInSecondsStop;
					float[] medianChromatogramArray;

					int correlationEncodedLength=rs.getInt(11);
					if (correlationEncodedLength==0) {
						// 0 indicates null, which indicates missing
						correlationArray=null;
						rtInSecondsStart=0.0f;
						rtInSecondsStop=0.0f;
						medianChromatogramArray=null;
					} else {
						correlationArray=ByteConverter.toFloatArray(CompressionUtils.decompress(rs.getBytes(12), correlationEncodedLength));
						rtInSecondsStart=rs.getFloat(13);
						rtInSecondsStop=rs.getFloat(14);
						int medianChromatogramEncodedLength=rs.getInt(15);
						medianChromatogramArray=ByteConverter.toFloatArray(CompressionUtils.decompress(rs.getBytes(16), medianChromatogramEncodedLength));
					}

					HashSet<String> accessions=PSMData.stringToAccessions(rs.getString(17));
					String sourceFile=rs.getString(18);
					if (correlationEncodedLength==0) {
						entry.add(new LibraryEntry(sourceFile, accessions, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray));
					} else {
						entry.add(new ChromatogramLibraryEntry(sourceFile, accessions, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray,
								correlationArray, medianChromatogramArray, new Range(rtInSecondsStart, rtInSecondsStop)));
					}
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
		Connection c=getConnection();
		try {
			Statement s=c.createStatement();
			try {
				try {
					Version version=getVersion();
					if (new Version(0, 1, 2).amIAbove(version)&&version.amIAbove(new Version(0, 0, 9))) {
						s.execute("ALTER TABLE entries ADD COLUMN CorrelationEncodedLength int");
						s.execute("ALTER TABLE entries ADD COLUMN CorrelationArray blob");
						s.execute("ALTER TABLE entries ADD COLUMN RTInSecondsStart double");
						s.execute("ALTER TABLE entries ADD COLUMN RTInSecondsStop double");
						s.execute("ALTER TABLE entries ADD COLUMN MedianChromatogramEncodedLength int");
						s.execute("ALTER TABLE entries ADD COLUMN MedianChromatogramArray blob");
					}
					if (new Version(0, 1, 4).amIAbove(version)&&version.amIAbove(new Version(0, 1, 2))) {
						s.execute("ALTER TABLE fragmentquants ADD COLUMN background double");
						s.execute("ALTER TABLE fragmentquants ADD COLUMN PeptideSeq string");
						s.execute("ALTER TABLE peptidequants ADD COLUMN PeptideSeq string");
					}

					if (new Version(0, 1, 5).amIAbove(version)&&version.amIAbove(new Version(0, 1, 2))) {
						s.execute("ALTER TABLE peptidequants ADD COLUMN LocalizationPeptideModSeq string");
						s.execute("ALTER TABLE peptidequants ADD COLUMN LocalizationScore double");
						s.execute("ALTER TABLE peptidequants ADD COLUMN NumberOfMods int");
						s.execute("ALTER TABLE peptidequants ADD COLUMN IsSiteSpecific boolean");
						s.execute("ALTER TABLE peptidequants ADD COLUMN RTInSecondsCenter double");
					}
				} catch (SQLException sqle) {
					// the metadata table is missing, so do nothing and create
					// it in the next line
				}

				s.execute("CREATE TABLE IF NOT EXISTS metadata ( "+"Key string not null, Value string not null, "+"PRIMARY KEY (Key) "+")");

				s.execute("CREATE TABLE IF NOT EXISTS entries ( "
						+"PrecursorMz double not null, PrecursorCharge int not null, PeptideModSeq string not null, PeptideSeq string not null, Copies int not null, RTInSeconds double not null, Score double not null, MassEncodedLength int not null, MassArray blob not null, IntensityEncodedLength int not null, IntensityArray blob not null, CorrelationEncodedLength int, CorrelationArray blob, RTInSecondsStart double, RTInSecondsStop double, MedianChromatogramEncodedLength int, MedianChromatogramArray blob, SourceFile string not null, "
						+"PRIMARY KEY (PrecursorCharge, PeptideModSeq, SourceFile), "+"FOREIGN KEY (PeptideSeq) REFERENCES proteins (PeptideSeq) "+")");

				s.execute("CREATE TABLE IF NOT EXISTS proteins ( "+"PeptideSeq string not null, ProteinAccessions string not null, "+"PRIMARY KEY (PeptideSeq) "+")");

				s.execute("CREATE TABLE IF NOT EXISTS peptidequants ( "
						+"PrecursorCharge int not null, PeptideModSeq string not null, PeptideSeq string not null, SourceFile string not null, LocalizationPeptideModSeq string, LocalizationScore double, NumberOfMods int, IsSiteSpecific boolean, RTInSecondsCenter double not null, RTInSecondsStart double not null, RTInSecondsStop double not null, TotalIntensity double not null, NumberOfQuantIons int not null, BestFragmentCorrelation double not null, BestFragmentDeltaMassPPM double not null, MedianChromatogramEncodedLength int not null, MedianChromatogramArray blob not null,"
						+"PRIMARY KEY (PrecursorCharge, PeptideModSeq, SourceFile), "
						+"FOREIGN KEY (PrecursorCharge, PeptideModSeq, SourceFile) REFERENCES entries (PrecursorCharge, PeptideModSeq, SourceFile) "+")");

				s.execute("CREATE TABLE IF NOT EXISTS fragmentquants ( "
						+"PrecursorCharge int not null, PeptideModSeq string not null, PeptideSeq string not null, SourceFile string not null, IonType string not null, FragmentMass double not null, Correlation double not null, Background double not null, DeltaMassPPM double not null, Intensity double not null, "
						+"FOREIGN KEY (PrecursorCharge, PeptideModSeq, SourceFile) REFERENCES entries (PrecursorCharge, PeptideModSeq, SourceFile) "+")");

				c.commit();
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	public void dropIndices() throws IOException, SQLException {
		Connection c=getConnection();
		try {
			Statement s=c.createStatement();
			try {

				s.execute("drop index if exists \"PeptideModSeq_Entries_index\"");
				s.execute("drop index if exists \"PeptideSeq_Entries_index\"");
				s.execute("drop index if exists \"PrecursorMz_Entries_index\"");
				s.execute("drop index if exists \"ProteinAccessions_Proteins_index\"");
				s.execute("drop index if exists \"PeptideModSeq_Peptides_index\"");
				s.execute("drop index if exists \"PeptideModSeq_Fragments_index\"");
				s.execute("drop index if exists \"PeptideSeq_Peptides_index\"");
				s.execute("drop index if exists \"PeptideSeq_Fragments_index\"");

				c.commit();
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	public void createIndices() throws IOException, SQLException {
		Connection c=getConnection();
		try {
			Statement s=c.createStatement();
			try {

				s.execute("create index if not exists \"PeptideModSeq_Entries_index\" on \"entries\" (\"PeptideModSeq\" ASC)");
				s.execute("create index if not exists \"PeptideSeq_Entries_index\" on \"entries\" (\"PeptideSeq\" ASC)");
				s.execute("create index if not exists \"PrecursorMz_Entries_index\" on \"entries\" (\"PrecursorMz\" ASC)");
				s.execute("create index if not exists \"ProteinAccessions_Proteins_index\" on \"proteins\" (\"ProteinAccessions\" ASC)");
				s.execute("create index if not exists \"PeptideModSeq_Peptides_index\" on \"peptidequants\" (\"PeptideModSeq\" ASC)");
				s.execute("create index if not exists \"PeptideModSeq_Fragments_index\" on \"fragmentquants\" (\"PeptideModSeq\" ASC)");
				s.execute("create index if not exists \"PeptideSeq_Peptides_index\" on \"peptidequants\" (\"PeptideSeq\" ASC)");
				s.execute("create index if not exists \"PeptideSeq_Fragments_index\" on \"fragmentquants\" (\"PeptideSeq\" ASC)");

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
