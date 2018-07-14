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
import java.util.HashSet;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideAccessionMatchingTrie;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.map.hash.TIntFloatHashMap;
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
	
	@SuppressWarnings("resource") // this is properly closed, Eclipse just can't follow the if/else logic
	public void getCopyEntriesToLibrary(LibraryFile library, Optional<TObjectFloatHashMap<String>> irtMap, File fastaFile, SearchParameters params) throws IOException, SQLException, DataFormatException {
		Logger.logLine("Reading BLIB file");
		
		ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select RefSpectraID, retentionTime from RetentionTimes where bestSpectrum=1");
				TIntFloatHashMap rtMap=new TIntFloatHashMap();
				while (rs.next()) {
					int refSpectraID=rs.getInt(1);
					float rt=rs.getFloat(2);
					rtMap.put(refSpectraID, rt);
				}

				boolean hasScore=doesColumnExist(tempFile, "RefSpectra", "score");
				boolean hasFileID=doesColumnExist(tempFile, "RefSpectra", "fileID");
				if (hasFileID) {
					if (hasScore) {
						rs=s.executeQuery(
								"select RefSpectra.id, RefSpectra.precursorMZ, RefSpectra.precursorCharge, RefSpectra.peptideModSeq, RefSpectra.copies, RefSpectra.numPeaks, RefSpectraPeaks.peakMZ, RefSpectraPeaks.peakIntensity, RefSpectra.score, SpectrumSourceFiles.fileName "
										+"from RefSpectra, RefSpectraPeaks, SpectrumSourceFiles "
										+"where RefSpectra.id == RefSpectraPeaks.RefSpectraID and SpectrumSourceFiles.id == RefSpectra.fileID");
					} else {
						rs=s.executeQuery(
								"select RefSpectra.id, RefSpectra.precursorMZ, RefSpectra.precursorCharge, RefSpectra.peptideModSeq, RefSpectra.copies, RefSpectra.numPeaks, RefSpectraPeaks.peakMZ, RefSpectraPeaks.peakIntensity, SpectrumSourceFiles.fileName "
										+"from RefSpectra, RefSpectraPeaks, SpectrumSourceFiles "
										+"where RefSpectra.id == RefSpectraPeaks.RefSpectraID and SpectrumSourceFiles.id == RefSpectra.fileID");
					}
				} else {
					if (hasScore) {
						rs=s.executeQuery(
								"select RefSpectra.id, RefSpectra.precursorMZ, RefSpectra.precursorCharge, RefSpectra.peptideModSeq, RefSpectra.copies, RefSpectra.numPeaks, RefSpectraPeaks.peakMZ, RefSpectraPeaks.peakIntensity, RefSpectra.score "
										+"from RefSpectra, RefSpectraPeaks "
										+"where RefSpectra.id == RefSpectraPeaks.RefSpectraID");
					} else {
						rs=s.executeQuery(
								"select RefSpectra.id, RefSpectra.precursorMZ, RefSpectra.precursorCharge, RefSpectra.peptideModSeq, RefSpectra.copies, RefSpectra.numPeaks, RefSpectraPeaks.peakMZ, RefSpectraPeaks.peakIntensity "
										+"from RefSpectra, RefSpectraPeaks, SpectrumSourceFiles "
										+"where RefSpectra.id == RefSpectraPeaks.RefSpectraID");
					}
				}
				
				AminoAcidConstants constants=new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());
				int missing=0;
				int total=0;
				while (rs.next()) {
					int refSpectraID=rs.getInt(1);
					double precursorMZ=rs.getDouble(2);
					byte precursorCharge=(byte)rs.getInt(3);
					String peptideModSeq=rs.getString(4);
					
					// precursors not set? This is a bug in Skyline exporting
					FragmentationModel model=new FragmentationModel(peptideModSeq, constants);
					precursorMZ=model.getChargedMass(precursorCharge);
					
					int copies=rs.getInt(5);
					int numPeaks=rs.getInt(6);
					double[] massArray=decompressDouble(rs.getBytes(7), numPeaks);
					float[] intensityArray=decompressFloat(rs.getBytes(8), numPeaks);
					float score;
					String sourceFile;
					if (hasFileID) {
						if (hasScore) {
							score=(float)rs.getDouble(9);
							sourceFile=rs.getString(10);
						} else {
							score=0.0f;
							sourceFile=rs.getString(9);
						}
					} else {
						if (hasScore) {
							score=(float)rs.getDouble(9);
							sourceFile="unknown";
						} else {
							score=0.0f;
							sourceFile="unknown";
						}
					}
					
					float retentionTime=rtMap.get(refSpectraID);
					if (irtMap.isPresent()) {
						if (irtMap.get().contains(peptideModSeq)) {
							retentionTime=irtMap.get().get(peptideModSeq);
						} else {
							missing++;
						}
					}
					// RT are usually in minutes, but not always depending on blib version. Warping 
					// won't be bothered by incorrect absolute values if this assumption is wrong.
					retentionTime=retentionTime*60.0f;
					total++;

					entries.add(new LibraryEntry(sourceFile, new HashSet<String>(), precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray));
				}
				if (missing>0) {
					Logger.logLine("Missing iRT for "+missing+" of "+total+" peptides, using RT in file.");
				}

				Logger.logLine("Reading Fasta file "+fastaFile.getName());
				ArrayList<FastaEntryInterface> proteins=FastaReader.readFasta(fastaFile);
				
				Logger.logLine("Constructing trie from library peptides");
				PeptideAccessionMatchingTrie trie=new PeptideAccessionMatchingTrie(entries);
				trie.addFasta(proteins);

				int[] counts=new int[21];
				for (LibraryEntry entry : entries) {
					int size=Math.min(counts.length-1, entry.getAccessions().size());
					counts[size]++;
				}
				Logger.logLine("Accession count histogram: ");
				for (int i=0; i<counts.length; i++) {
					Logger.logLine(i+" Acc\t"+counts[i]+" Counts");
				}

				if (counts[0]>0) {
					Logger.errorLine(counts[0]+" library entries can't be linked to proteins! These entries will be dropped.");
				}
				Logger.logLine("Writing library file "+library.getName());
				library.dropIndices();
				library.addEntries(entries);
				library.createIndices();
				
				rs.close();
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
		
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
		TCharDoubleHashMap fixedMods=aaConstants.getFixedMods();
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
			PreparedStatement prepRTs=c.prepareStatement("insert into RetentionTimes (RefSpectraID, RedundantRefSpectraID, SpectrumSourceID, retentionTime, startTime, endTime, bestSpectrum) VALUES (?,?,?,?,?,?,?)");
			PreparedStatement prepMods=c.prepareStatement("insert into Modifications (id, RefSpectraID, position, mass) VALUES (?,?,?,?)");
			
			try {
				for (LibraryEntry entry : entries) {
					idCounter++;
					
					double[] massArray;
					float[] intensityArray;
					Range rtRange;
					if (entry instanceof IntegratedLibraryEntry) {
						TransitionRefinementData data=((IntegratedLibraryEntry) entry).getRefinementData();
						rtRange=((IntegratedLibraryEntry)entry).getRtRange();
						Optional<double[]> masses=data.getMassArray();
						if (masses.isPresent()) {
							massArray=masses.get();
							intensityArray=data.getIntensityArray().get();
						} else {
							massArray=entry.getMassArray();
							intensityArray=entry.getIntensityArray();
						}
					} else {
						rtRange=new Range(entry.getRetentionTime()-60, entry.getRetentionTime()+60);
						massArray=entry.getMassArray();
						intensityArray=entry.getIntensityArray();
					}
					
					prep.setInt(1, idCounter); // id
					prep.setString(2, entry.getPeptideSeq()); // pepSeq
					prep.setDouble(3, entry.getPrecursorMZ()); // precursorMZ
					prep.setInt(4, entry.getPrecursorCharge()); // precursorCharge
					String peptideModSeq=PeptideUtils.formatForSkyline(entry.getPeptideModSeq(), aaConstants);
					prep.setString(5, peptideModSeq); // peptideModSeq
					prep.setString(6, "-"); // prevAA
					prep.setString(7, "-"); // nextAA
					prep.setInt(8, 1); // copies
					prep.setInt(9, massArray.length); // numPeaks
					prep.setDouble(10, entry.getRetentionTime()/60f); // retentionTime
					prep.setInt(11, jobCounter); // fileID
					prep.setString(12, diaFileName+"."+entry.getSpectrumIndex()+"."+entry.getSpectrumIndex()+"."+entry.getPrecursorCharge()); // SpecIDinFile
					prep.setDouble(13, entry.getScore()); // score
					prep.setByte(14, scoreTypeID); // scoreType
					prep.addBatch();
					
					prepPeaks.setInt(1,  idCounter);
					prepPeaks.setBytes(2, compressDouble(massArray));
					prepPeaks.setBytes(3, compressFloat(intensityArray));
					prepPeaks.addBatch();
					
					prepRTs.setInt(1,  idCounter);
					prepRTs.setInt(2,  0);
					prepRTs.setInt(3,  jobCounter);
					prepRTs.setDouble(4,  entry.getRetentionTime()/60f); // convert to minutes
					prepRTs.setDouble(5, rtRange.getStart()/60f);
					prepRTs.setDouble(6, rtRange.getStop()/60f);
					prepRTs.setInt(7,  1);
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
							prepMods.setDouble(4, fixedMods.get(aas[i].charAt(0)));
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
							prepMods.setDouble(4, mass);
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
		byte[] uncompressed=ByteConverter.toByteArray(masses, ByteOrder.LITTLE_ENDIAN);
		byte[] compressed=CompressionUtils.compress(uncompressed);
		if (compressed.length<uncompressed.length) {
			return compressed;
		} else {
			return uncompressed;
		}
	}
	
	private byte[] compressFloat(float[] intensities) throws IOException {
		byte[] uncompressed=ByteConverter.toByteArray(intensities, ByteOrder.LITTLE_ENDIAN);
		byte[] compressed=CompressionUtils.compress(uncompressed);
		if (compressed.length<uncompressed.length) {
			return compressed;
		} else {
			return uncompressed;
		}
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
				s.execute("CREATE TABLE if not exists RetentionTimes (RefSpectraID INTEGER, RedundantRefSpectraID INTEGER, SpectrumSourceID INTEGER, ionMobilityValue REAL, ionMobilityType INTEGER, retentionTime REAL, startTime REAL, endTime REAL, bestSpectrum INTEGER, FOREIGN KEY(RefSpectraID) REFERENCES RefSpectra(id)) ");
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
