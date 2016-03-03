package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import gnu.trove.map.hash.TCharFloatHashMap;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class BlibFile extends SQLFile {
	public static final String BLIB=".blib";
	private final File tempFile;
	private File userFile;

	public BlibFile() throws IOException {
		tempFile=File.createTempFile("encyclopedia_", BLIB);
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
	
	public void setUserFile(File userFile) {
		this.userFile=userFile;
	}
	
	public void getStreamEntriesToLibrary(LibraryFile library, Optional<TObjectFloatHashMap<String>> irtMap) throws IOException, SQLException, DataFormatException {
		library.dropIndices();
		
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s
						.executeQuery("select RefSpectra.precursorMZ, RefSpectra.precursorCharge, RefSpectra.peptideModSeq, RefSpectra.copies, RefSpectra.numPeaks, RefSpectra.retentionTime, RefSpectra.score, RefSpectraPeaks.peakMZ, RefSpectraPeaks.peakIntensity from RefSpectra, RefSpectraPeaks "
								+"where RefSpectra.id == RefSpectraPeaks.RefSpectraID");

				ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
				int missing=0;
				int total=0;
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
					
					if (irtMap.isPresent()) {
						if (irtMap.get().contains(peptideModSeq)) {
							retentionTime=irtMap.get().get(peptideModSeq);
						} else {
							missing++;
						}
					}
					total++;

					entries.add(new LibraryEntry(precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray));
					
					if (entries.size()>1000) {
						library.addEntries(entries);
						entries.clear();
						Logger.log(".");
					}
				}
				if (missing>0) {
					System.out.println("Missing iRT for "+missing+" of "+total+" peptides, using RT in file.");
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
	public int[] addLibrary(SearchJobData job, ArrayList<LibraryEntry> entries, int idCounter, int jobCounter, int modCounter) throws IOException, SQLException {
		String diaFileName=job.getDiaFile().getName();
		String spectrumIDPrefix=diaFileName;
		if (spectrumIDPrefix.indexOf('.')>0) {
			spectrumIDPrefix=spectrumIDPrefix.substring(0, spectrumIDPrefix.indexOf('.'));
		}
		
		String rootName=userFile.getName();
		if (rootName.endsWith(BLIB)) {
			rootName=rootName.substring(0, rootName.length()-5);
		}
		
		AminoAcidConstants aaConstants=job.getParameters().getAAConstants();
		TCharFloatHashMap fixedMods=aaConstants.getFixedMods();
		char[] fixedModdedAAs=fixedMods.keys();
		modCounter++;
		
		Connection c=getConnection(tempFile);
		try {
			Statement normalStatement=c.createStatement();
			jobCounter++;
			normalStatement.executeUpdate("insert into SpectrumSourceFiles (id, fileName) VALUES ("+jobCounter+",\""+diaFileName+"\");");
			
			String libLSID="urn:lsid:proteome.gs.washington.edu:spectral_library:pecan:"+rootName;
			SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
			int totalSpectra=entries.size()+idCounter;
			normalStatement.executeUpdate("insert into LibInfo (libLSID, createTime, numSpecs, majorVersion, minorVersion) VALUES ("+
					"\""+libLSID+"\",\""+format.format(new Date())+"\","+totalSpectra+",1,"+jobCounter+");");
			
			ResultSet results=normalStatement.executeQuery("select count(*) from ScoreTypes");
			int numberOfScores=results.getInt(1);
			results.close();
			
			byte scoreTypeID=1;
			if (numberOfScores==0) {
				normalStatement.executeUpdate("insert into ScoreTypes (id, scoreType) VALUES ("+scoreTypeID+",\"Pecan_"+job.getVersion()+"\");");
			}
			
			normalStatement.close();
			
			PreparedStatement prep=c.prepareStatement("insert into RefSpectra (id, peptideSeq, precursorMZ, precursorCharge, peptideModSeq, prevAA, nextAA, copies, numPeaks, retentionTime, fileID, SpecIDinFile, score, scoreType) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
			PreparedStatement prepPeaks=c.prepareStatement("insert into RefSpectraPeaks (RefSpectraID, peakMZ, peakIntensity) VALUES (?,?,?)");
			PreparedStatement prepRTs=c.prepareStatement("insert into RetentionTimes (RefSpectraID, RedundantRefSpectraID, SpectrumSourceID, retentionTime, bestSpectrum) VALUES (?,?,?,?,?)");
			PreparedStatement prepMods=c.prepareStatement("insert into Modifications (id, RefSpectraID, position, mass) VALUES (?,?,?,?)");
			
			try {
				for (LibraryEntry entry : entries) {
					idCounter++;
					prep.setInt(1, idCounter); // id
					prep.setString(2, entry.getPeptideSeq()); // pepSeq
					prep.setDouble(3, entry.getPrecursorMZ()); // precursorMZ
					prep.setInt(4, entry.getPrecursorCharge()); // precursorCharge
					String peptideModSeq=entry.getPeptideModSeq();
					for (char aa : fixedModdedAAs) {
						float mass=fixedMods.get(aa);
						String replacement=aa+(mass>=0?"[+":"[")+mass+"]";
						peptideModSeq=peptideModSeq.replace(Character.toString(aa), replacement);
					}
					prep.setString(5, peptideModSeq); // peptideModSeq
					prep.setString(6, "-"); // prevAA
					prep.setString(7, "-"); // nextAA
					prep.setInt(8, 1); // copies
					prep.setInt(9, entry.getMassArray().length); // numPeaks
					prep.setDouble(10, entry.getRetentionTime()/60f); // retentionTime
					prep.setInt(11, jobCounter); // fileID
					prep.setString(12, diaFileName+"."+entry.getSpectrumIndex()+"."+entry.getSpectrumIndex()+"."+entry.getPrecursorCharge()); // SpecIDinFile
					prep.setDouble(13, entry.getScore()); // score
					prep.setByte(14, scoreTypeID); // scoreType
					prep.addBatch();
					
					prepPeaks.setInt(1,  idCounter);
					prepPeaks.setBytes(2, compressDouble(entry.getMassArray()));
					prepPeaks.setBytes(3, compressFloat(entry.getIntensityArray()));
					prepPeaks.addBatch();
					
					prepRTs.setInt(1,  idCounter);
					prepRTs.setInt(2,  0);
					prepRTs.setInt(3,  jobCounter);
					prepRTs.setDouble(4,  entry.getRetentionTime()/60f); // convert to minutes
					prepRTs.setInt(5,  0);
					prepRTs.addBatch();
					
					FragmentationModel model=new FragmentationModel(peptideModSeq, aaConstants);
					String[] aas=model.getAas();
					for (int i=0; i<aas.length; i++) {
						boolean added=false;
						Pair<Character, Double> aa=FragmentationModel.parseAA(aas[i]);
						if (fixedMods.contains(aa.x)) {
							prepMods.setInt(1, modCounter);
							prepMods.setInt(2, idCounter);
							prepMods.setInt(3, (i+1));
							prepMods.setFloat(4, fixedMods.get(aas[i].charAt(0)));
							added=true;
						} else if (aa.y!=null) {
							float mass=aa.y.floatValue();
							int index;
							if (aas[i].charAt(0)=='[') {
								// prefix mod
								index=i;
							} else {
								// post mod
								index=i+1;
							}
							prepMods.setInt(1, modCounter);
							prepMods.setInt(2, idCounter);
							prepMods.setInt(3, index);
							prepMods.setFloat(4, mass);
							added=true;
						}
						if (added) {
							modCounter++;
							prepMods.addBatch();
						}
					}
				}
				prep.executeBatch();
				prep.close();

				prepPeaks.executeBatch();
				prepPeaks.close();

				prepRTs.executeBatch();
				prepRTs.close();

				prepMods.executeBatch();
				prepMods.close();
				
				c.commit();
				
				return new int[] {idCounter, jobCounter, modCounter};
				
			} finally {
				prep.close();
			}
		} finally {
			c.close();
		}
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
	
	private byte[] compressDouble(double[] masses) throws IOException {
		return CompressionUtils.compress(ByteConverter.toByteArray(masses, ByteOrder.LITTLE_ENDIAN));
	}
	
	private byte[] compressFloat(float[] intensities) throws IOException {
		return CompressionUtils.compress(ByteConverter.toByteArray(intensities, ByteOrder.LITTLE_ENDIAN));
	}
	
	private void createNewTables() throws IOException, SQLException {
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				s.execute("CREATE TABLE if not exists LibInfo(libLSID TEXT, createTime TEXT, numSpecs INTEGER, majorVersion INTEGER, minorVersion INTEGER)");
				s.execute("CREATE TABLE if not exists Modifications (id INTEGER primary key autoincrement not null,RefSpectraID INTEGER, position INTEGER, mass REAL)");
				s.execute("CREATE TABLE if not exists RefSpectra (id INTEGER primary key autoincrement not null, peptideSeq VARCHAR(150), precursorMZ REAL, precursorCharge INTEGER, peptideModSeq VARCHAR(200), prevAA CHAR(1), nextAA CHAR(1), copies INTEGER, numPeaks INTEGER, ionMobilityValue REAL, ionMobilityType INTEGER, retentionTime REAL, fileID INTEGER, SpecIDinFile VARCHAR(256), score REAL, scoreType TINYINT)");
				s.execute("CREATE TABLE if not exists RefSpectraPeaks(RefSpectraID INTEGER, peakMZ BLOB, peakIntensity BLOB)");
				s.execute("CREATE TABLE if not exists RetentionTimes (RefSpectraID INTEGER, RedundantRefSpectraID INTEGER, SpectrumSourceID INTEGER, ionMobilityValue REAL, ionMobilityType INTEGER, retentionTime REAL, bestSpectrum INTEGER, FOREIGN KEY(RefSpectraID) REFERENCES RefSpectra(id)) ");
				s.execute("CREATE TABLE if not exists ScoreTypes (id INTEGER PRIMARY KEY, scoreType VARCHAR(128)) ");
				s.execute("CREATE TABLE if not exists SpectrumSourceFiles (id INTEGER PRIMARY KEY autoincrement not null,fileName VARCHAR(512)) ");
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
}
