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
import java.util.zip.DataFormatException;


import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Chromatogram;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ChromatogramLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.Version;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class LibraryFile extends SQLFile implements LibraryInterface {
	public static final String ELIB=".elib";
	public static final String VERSION_STRING="version";
	public static final Version[] ACCEPTABLE_VERSIONS=new Version[] {new Version(0, 1, 0), new Version(0, 1, 1), new Version(0, 1, 2)};
	public static final Version MOST_RECENT_VERSION=new Version(0, 1, 2);
	
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
		map.put(VERSION_STRING, MOST_RECENT_VERSION.toString());
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
	
	public Version getVersion() throws IOException, SQLException {
		HashMap<String, String> meta=getMetadata();
		return new Version(meta.get(VERSION_STRING));
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
		Connection c=getConnection(tempFile);
		try {
			PreparedStatement peptidePrep=c.prepareStatement("INSERT INTO peptidequants (PrecursorCharge, PeptideModSeq, SourceFile, RTInSecondsStart, RTInSecondsStop, TotalIntensity, NumberOfQuantIons, BestFragmentCorrelation, BestFragmentDeltaMassPPM, MedianChromatogramEncodedLength, MedianChromatogramArray) VALUES (?,?,?,?,?,?,?,?,?,?,?)");
			PreparedStatement fragmentPrep=c.prepareStatement("INSERT INTO fragmentquants (PrecursorCharge, PeptideModSeq, SourceFile, IonType, FragmentMass, Correlation, DeltaMassPPM, Intensity) VALUES (?,?,?,?,?,?,?,?)");
			try {
				for (LibraryEntry recast : uniqueEntries) {
					IntegratedLibraryEntry entry=(IntegratedLibraryEntry)recast;
					TransitionRefinementData data=entry.getRefinementData();

					float[] correlationArray=data.getCorrelationArray();
					float[] integrationArray=data.getIntegrationArray();
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
					
					peptidePrep.setInt(1, entry.getPrecursorCharge());
					peptidePrep.setString(2, entry.getPeptideModSeq());
					peptidePrep.setString(3, entry.getSource());
					peptidePrep.setFloat(4, data.getRange().getStart());
					peptidePrep.setFloat(5, data.getRange().getStop());
					peptidePrep.setFloat(6, data.getTotalIntensity(minimumCorrelation));
					peptidePrep.setInt(7, data.getTotalQuantIons(minimumCorrelation));
					peptidePrep.setFloat(8, bestCorrelation);
					peptidePrep.setFloat(9, bestDeltaMass);
					byte[] intensityByteArray=ByteConverter.toByteArray(data.getMedianChromatogram());
					peptidePrep.setInt(10, intensityByteArray.length);
					peptidePrep.setBytes(11, CompressionUtils.compress(intensityByteArray));
					peptidePrep.addBatch();
					
					for (int i=0; i<correlationArray.length; i++) {
						if (correlationArray[i]>=minimumCorrelation) {
							fragmentPrep.setInt(1, entry.getPrecursorCharge());
							fragmentPrep.setString(2, entry.getPeptideModSeq());
							fragmentPrep.setString(3, entry.getSource());
							fragmentPrep.setString(4, "BY");
							fragmentPrep.setDouble(5, fragmentMassArray[i]);
							fragmentPrep.setFloat(6, correlationArray[i]);
							fragmentPrep.setFloat(7, ppmArray[i]);
							fragmentPrep.setFloat(8, integrationArray[i]);
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
			PreparedStatement prep=c.prepareStatement("INSERT INTO entries (PrecursorMZ, PrecursorCharge, PeptideModSeq, PeptideSeq, Copies, RTInSeconds, Score, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray, CorrelationEncodedLength, CorrelationArray, RTInSecondsStart, RTInSecondsStop, MedianChromatogramEncodedLength, MedianChromatogramArray, SourceFile) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
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
					
					
					if (entry instanceof Chromatogram) {
						Chromatogram cast=(Chromatogram)entry;
						
						byte[] correlationByteArray=ByteConverter.toByteArray(cast.getCorrelationArray());
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

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface#getEntries(edu.washington.gs.maccoss.encyclopedia.datastructures.Range)
	 */
	@Override
	public ArrayList<LibraryEntry> getEntries(String peptideModSeq, byte charge, boolean sqrt) throws IOException, SQLException, DataFormatException {
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select e.PrecursorMZ, e.PrecursorCharge, e.PeptideModSeq, e.Copies, e.RTInSeconds, e.Score, e.MassEncodedLength, e.MassArray, e.IntensityEncodedLength, e.IntensityArray, e.CorrelationEncodedLength, e.CorrelationArray blob, e.RTInSecondsStart, e.RTInSecondsStop, e.MedianChromatogramEncodedLength, e.MedianChromatogramArray, p.ProteinAccessions, e.SourceFile from entries e, proteins p"
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
				ResultSet rs=s.executeQuery("select e.PrecursorMZ, e.PrecursorCharge, e.PeptideModSeq, e.Copies, e.RTInSeconds, e.Score, e.MassEncodedLength, e.MassArray, e.IntensityEncodedLength, e.IntensityArray, e.CorrelationEncodedLength, e.CorrelationArray blob, e.RTInSecondsStart, e.RTInSecondsStop, e.MedianChromatogramEncodedLength, e.MedianChromatogramArray, p.ProteinAccessions, e.SourceFile from entries e, proteins p"
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
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select e.PrecursorMZ, e.PrecursorCharge, e.PeptideModSeq, e.Copies, e.RTInSeconds, e.Score, e.MassEncodedLength, e.MassArray, e.IntensityEncodedLength, e.IntensityArray, e.CorrelationEncodedLength, e.CorrelationArray blob, e.RTInSecondsStart, e.RTInSecondsStop, e.MedianChromatogramEncodedLength, e.MedianChromatogramArray, p.ProteinAccessions, e.SourceFile from entries e, proteins p"
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
		Connection c=getConnection(tempFile);
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
				} catch (SQLException sqle) {
					// the metadata table is missing, so do nothing and create it in the next line
				}
				
				s.execute("CREATE TABLE IF NOT EXISTS metadata ( "
						+ "Key string not null, Value string not null, "
						+ "PRIMARY KEY (Key) "
						+ ")");
				
				s.execute("CREATE TABLE IF NOT EXISTS entries ( "
						+ "PrecursorMz double not null, PrecursorCharge int not null, PeptideModSeq string not null, PeptideSeq string not null, Copies int not null, RTInSeconds double not null, Score double not null, MassEncodedLength int not null, MassArray blob not null, IntensityEncodedLength int not null, IntensityArray blob not null, CorrelationEncodedLength int, CorrelationArray blob, RTInSecondsStart double, RTInSecondsStop double, MedianChromatogramEncodedLength int, MedianChromatogramArray blob, SourceFile string not null, "
						+ "PRIMARY KEY (PrecursorCharge, PeptideModSeq, SourceFile), "
						+ "FOREIGN KEY (PeptideSeq) REFERENCES proteins (PeptideSeq) "
						+ ")");
				
				s.execute("CREATE TABLE IF NOT EXISTS proteins ( "
						+ "PeptideSeq string not null, ProteinAccessions string not null, "
						+ "PRIMARY KEY (PeptideSeq) "
						+ ")");
				
				s.execute("CREATE TABLE IF NOT EXISTS peptidequants ( "
						+ "PrecursorCharge int not null, PeptideModSeq string not null, SourceFile string not null, RTInSecondsStart double not null, RTInSecondsStop double not null, TotalIntensity double not null, NumberOfQuantIons int not null, BestFragmentCorrelation double not null, BestFragmentDeltaMassPPM double not null, MedianChromatogramEncodedLength int not null, MedianChromatogramArray blob not null,"
						+ "PRIMARY KEY (PrecursorCharge, PeptideModSeq, SourceFile), "
						+ "FOREIGN KEY (PrecursorCharge, PeptideModSeq, SourceFile) REFERENCES entries (PrecursorCharge, PeptideModSeq, SourceFile) "
						+ ")");
				
				s.execute("CREATE TABLE IF NOT EXISTS fragmentquants ( "
						+ "PrecursorCharge int not null, PeptideModSeq string not null, SourceFile string not null, IonType string not null, FragmentMass double not null, Correlation double not null, DeltaMassPPM double not null, Intensity double not null, "
						+ "FOREIGN KEY (PrecursorCharge, PeptideModSeq, SourceFile) REFERENCES entries (PrecursorCharge, PeptideModSeq, SourceFile) "
						+ ")");
				
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
