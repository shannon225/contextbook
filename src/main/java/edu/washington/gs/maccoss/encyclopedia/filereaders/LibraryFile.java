package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.sql.Types;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.StringTokenizer;
import java.util.function.Predicate;
import java.util.zip.DataFormatException;

import com.google.common.collect.ImmutableList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeakLocationInferrerInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorProteinGroup;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.AmbiguousPeptideModSeq;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Chromatogram;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ChromatogramLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LocalizedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroupInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SimplePeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.io.Version;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Ion;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.QuantitativeDIAData;
import edu.washington.gs.maccoss.encyclopedia.utils.math.BenjaminiHochberg;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;

public class LibraryFile extends SQLFile implements LibraryInterface {
	public static boolean OPEN_IN_PLACE=false;

	private static final String SOURCEFILE_TIC_PREFIX="TIC_";
	private static final String SOURCEFILE_STRING="sourcefile";
	private static final String SOURCE_FILE_SPLIT="|";
	public static final String DLIB=".dlib";
	public static final String ELIB=".elib";
	public static final String ENCYCLOPEDIA_VERSION = "EncyclopediaVersion";
	public static final String PERCOLATOR_VERSION = "PercolatorVersion";
	public static final String UNKNOWN = "Unknown";
	public static final Version[] ACCEPTABLE_VERSIONS = new Version[] {
			new Version(0, 1, 0), new Version(0, 1, 1), new Version(0, 1, 2), new Version(0, 1, 3), new Version(0, 1, 4),
			new Version(0, 1, 5), new Version(0, 1, 6), new Version(0, 1, 7), new Version(0, 1, 8), new Version(0, 1, 9),
			new Version(0, 1, 10), new Version(0, 1, 11), new Version(0, 1, 12), new Version(0, 1, 13), new Version(0, 1, 14)
	};
	public static final Version MOST_RECENT_VERSION=new Version(0, 1, 14);

	private File userFile=null;
	private File tempFile;

	public LibraryFile() throws IOException {
		tempFile=File.createTempFile("encyclopedia_", ELIB);
		tempFile.deleteOnExit();
	}

	public static boolean isVersionAcceptable(Version version) {
		if (version==null)
			return false;

		for (Version string : ACCEPTABLE_VERSIONS) {
			if (string.equals(version))
				return true;
		}
		return false;
	}

	public String getName() {
		return userFile==null ? tempFile.getName() : userFile.getName();
	}

	public void openFile(File userFile) throws IOException, SQLException {
		this.userFile=userFile;
		openFile();
	}

	public void openFile() throws IOException, SQLException {
		if (OPEN_IN_PLACE) {
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
		String version = getClass().getPackage().getImplementationVersion();
		addMetadata(ENCYCLOPEDIA_VERSION, version == null ? UNKNOWN : version);
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
		String key=SOURCEFILE_TIC_PREFIX+ diaFile.getOriginalFileName();

		HashMap<String, String> map=new HashMap<String, String>();
		map.put(key, Float.toString(diaFile.getTIC()));

		addMetadata(map);
	}

	public float getTIC(StripeFileInterface diaFile) throws IOException, SQLException {
		return getTIC(diaFile.getOriginalFileName());
	}

	public float getTIC(String originalFileName) throws IOException, SQLException {
		String key=SOURCEFILE_TIC_PREFIX+originalFileName;

		String value=getMetadata().get(key);
		if (value==null)
			return 0.0f;
		return Float.parseFloat(value);
	}

	public void addRtAlignment(SearchJobData job, PeakLocationInferrerInterface inferrer) {
		Optional.ofNullable(inferrer.getAlignmentData(job))
				.ifPresent(alignment -> addRtAlignment(job, alignment));
	}

	public void addRtAlignment(SearchJobData job, List<RetentionTimeAlignmentInterface.AlignmentDataPoint> alignment) {
		final String sourceFile = job.getDiaFileReader().getOriginalFileName();

		try (Connection c = getConnection()) {
			c.setAutoCommit(false);

			try (PreparedStatement s = c.prepareStatement(
					"INSERT INTO retentiontimes (SourceFile, Library, Actual, Predicted, Delta, Probability, Decoy, PeptideModSeq) VALUES (?, ?, ?, ?, ?, ?, ?, ?)"
			)) {
				for (int i = 0, alignmentSize = alignment.size(); i < alignmentSize; i++) {
					final RetentionTimeAlignmentInterface.AlignmentDataPoint pt = alignment.get(i);

					s.setString(1, sourceFile);

					// Convert values in minutes to seconds
					s.setFloat(2, 60*pt.getLibrary());
					s.setFloat(3, 60*pt.getActual());
					s.setFloat(4, 60*pt.getPredictedActual());
					s.setFloat(5, 60*pt.getDelta());

					s.setFloat(6, pt.getProbability());
					s.setObject(7, pt.isDecoy()); // used setObject to handle null Boolean
					s.setString(8, pt.getPeptideModSeq());

					s.addBatch();

					if (i % 8192 == 0) {
						s.executeBatch();
						s.clearBatch();
					}
				}
				s.executeBatch();
			} catch (SQLException e) {
				c.rollback();
			} finally {
				c.commit();
			}
		} catch (SQLException | IOException e) {
			e.printStackTrace();
		}
	}

	public void setSources(List<? extends SearchJobData> sources) throws IOException, SQLException {
		HashMap<String, String> map=new HashMap<String, String>();
		StringBuilder sb=new StringBuilder();
		for (SearchJobData searchJobData : sources) {
			if (sb.length()>0) {
				sb.append(SOURCE_FILE_SPLIT);
			}
			
			sb.append(searchJobData.getDiaFileReader().getOriginalFileName());

		}
		map.put(SOURCEFILE_STRING, sb.toString());
		addMetadata(map);
	}

	public void setFileVersion() throws IOException, SQLException {
		HashMap<String, String> map=new HashMap<String, String>();
		map.put(VERSION_STRING, MOST_RECENT_VERSION.toString());
		addMetadata(map);
	}

	public void addMetadata(String key, String value) throws IOException, SQLException {
		HashMap<String, String> data=new HashMap<>();
		data.put(key, value);
		addMetadata(data);
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

	public List<Path> getSourceFiles() throws IOException, SQLException {
		final String sources=getMetadata().get(SOURCEFILE_STRING);

		if (null == sources) {
			return Collections.emptyList();
		}

		final StringTokenizer st=new StringTokenizer(sources, SOURCE_FILE_SPLIT);

		ImmutableList.Builder<Path> builder = ImmutableList.builder();
		while (st.hasMoreTokens()) {
			builder.add(Paths.get(st.nextToken()));
		}
		return builder.build();
	}

	public Optional<Path> getSource(SearchParameters parameters) {
		try {
			final Collection<Path> sources=getSourceFiles();

			if (sources.size()==0||sources.size()>1) {
				return Optional.empty();
			}

			Path next=sources.iterator().next();
			File file=next.toFile();
			if (file!=null) {
				if (file.exists()) {
					return Optional.ofNullable(next);
				} else if (new File(getFile().getParent(), file.getName()).exists()) {
					return Optional.ofNullable(new File(getFile().getParent(), file.getName()).toPath());
				}
			}
			return Optional.ofNullable(next);
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

	/**
	 * if there are repeated integrated peptides, we prefer peptides with higher scores
	 * @param entries
	 * @param inferrer
	 * @param localizationData
	 * @param aaConstants
	 * @param fdrThreshold
	 * @throws IOException
	 * @throws SQLException
	 */
	public void addIntegratedEntries(ArrayList<IntegratedLibraryEntry> entries, Optional<PeakLocationInferrerInterface> inferrer, Optional<HashMap<String, ModificationLocalizationData>> localizationData, AminoAcidConstants aaConstants, float fdrThreshold)
			throws IOException, SQLException {
		// first add normal data
		HashMap<String, LibraryEntry> repeatsCatcher=new HashMap<String, LibraryEntry>();
		HashMap<String, String> ptmRepeatsCatcher=new HashMap<String, String>();

        TObjectDoubleHashMap<String> localizationFDRValues=new TObjectDoubleHashMap<>();
        if (localizationData.isPresent()) {
            ArrayList<String> peptideModSeqs=new ArrayList<>();
            TDoubleArrayList pValues=new TDoubleArrayList();
            for (Entry<String, ModificationLocalizationData> entry : localizationData.get().entrySet()) {
            		peptideModSeqs.add(entry.getKey());
                double pvalue=Math.pow(10, -entry.getValue().getLocalizationScore());
                pValues.add(pvalue);
            }
            double[] fdrValues=BenjaminiHochberg.calculateAdjustedPValues(pValues.toArray());
            for (int i=0; i<fdrValues.length; i++) {
                localizationFDRValues.put(peptideModSeqs.get(i), fdrValues[i]);
            }
        }

		int countSkipped=0;
		for (IntegratedLibraryEntry entry : entries) {
			if (localizationData.isPresent()) {
				if(!localizationFDRValues.contains(entry.getPeptideModSeq())) {
					System.out.println("Skipping "+entry.getPeptideModSeq()+", +not localized!");
					countSkipped++;
					continue;
				}
				double fdr=localizationFDRValues.get(entry.getPeptideModSeq());
				if (fdr>fdrThreshold) {
					System.out.println("Skipping "+entry.getPeptideModSeq()+", FDR over threshold:"+fdr+"!");
					countSkipped++;
					continue;
				}	
			}
			String key=entry.getPeptideModSeq()+"+"+entry.getPrecursorCharge()+","+entry.getSource();
			LibraryEntry prev=repeatsCatcher.get(key);
			if (prev==null) {
				repeatsCatcher.put(key, entry);
			} else {
				// Logger.errorLine("Found collision writing elib: "+key+"
				// ("+entry.getScore()+" vs"+prev.getScore()+"), keeping best
				// scoring. Let Brian know if you see this!");
				if (entry.getScore()>prev.getScore()) {
					// prefer peptides with higher scores
					repeatsCatcher.put(key, entry);
				}
			}
		}
		if (countSkipped>0) {
			Logger.logLine("Skipped "+countSkipped+" peptides that did not pass localization FDR");
		}
		ArrayList<LibraryEntry> uniqueEntries=new ArrayList<LibraryEntry>(repeatsCatcher.values());
		Logger.logLine("Writing "+uniqueEntries.size()+" peptides to entries table...");
		addEntries(uniqueEntries);

		// then add integrated areas
		Connection c=getConnection();
		try {
			ArrayList<Pair<TransitionRefinementData, String>> dataAndSourceList=new ArrayList<Pair<TransitionRefinementData, String>>();
			for (LibraryEntry recast : uniqueEntries) {
				IntegratedLibraryEntry entry=(IntegratedLibraryEntry) recast;
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
							int prevAmbiguityScore=prevLocalizationData.isPresent() ? prevLocalizationData.get().getLocalizationPeptideModSeq().getNumAmbigousResidues() : 0;
							int newAmbiguityScore=value.getLocalizationData().isPresent() ? value.getLocalizationData().get().getLocalizationPeptideModSeq().getNumAmbigousResidues() : 0;
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

				for (Entry<String, TransitionRefinementData> mapEntry : uniqueDataMap.entrySet()) {
					TransitionRefinementData uniqueData=mapEntry.getValue();

					String key=uniqueData.getPeptideModSeq()+"+"+uniqueData.getPrecursorCharge()+","+sourceFile;
					if (ptmRepeatsCatcher.containsKey(key)) {
						System.err.println("FOUND EXTERNAL COLLISION, SKIPPING! "+uniqueData.getPeptideModSeq()+" (from:"+mapEntry.getKey()+")");
						System.err.println("PREV: "+ptmRepeatsCatcher.get(key));
						if (uniqueData.getLocalizationData().isPresent()) {
							ModificationLocalizationData modData=uniqueData.getLocalizationData().get();
							System.err.println("NEW:  "+modData.getLocalizationPeptideModSeq().getPeptideAnnotation());
						} else {
							System.err.println("NEW:  NO LOC: "+uniqueData.getPeptideModSeq());
						}
						continue;

					} else {
						if (uniqueData.getLocalizationData().isPresent()) {
							ModificationLocalizationData modData=uniqueData.getLocalizationData().get();
							ptmRepeatsCatcher.put(key, modData.getLocalizationPeptideModSeq().getPeptideAnnotation());
						} else {
							ptmRepeatsCatcher.put(key, "NO LOC: "+uniqueData.getPeptideModSeq());
						}
					}
					dataAndSourceList.add(new Pair<TransitionRefinementData, String>(uniqueData, sourceFile));
				}

			}

			Logger.logLine("Writing "+dataAndSourceList.size()+" peptides to peptidequants table...");

			// Issue 25 - skip entries that the RT inferrer could not process.
			// TODO: should this throw if too many entries are skipped?
			// TODO: instead of skipping these entries, should something else be
			// done?
			int numFilteredEntries=0;
			if (inferrer.isPresent()) {
				ArrayList<Pair<TransitionRefinementData, String>> filteredDataAndSourceList=new ArrayList<Pair<TransitionRefinementData, String>>();

				for (Pair<TransitionRefinementData, String> pair : dataAndSourceList) {
					if (inferrer.get().getTopNBestIons(pair.x.getPeptideModSeq(), pair.x.getPrecursorCharge())!=null) {
						filteredDataAndSourceList.add(pair);
					} else {
						numFilteredEntries++;
					}

				}

				dataAndSourceList=filteredDataAndSourceList;

				if (numFilteredEntries>0) {
					Logger.errorLine("Skipped "+numFilteredEntries+" integrated library entries because the RT inferrer could not find any fragment ions.");
				}
			}

			int start=0;
			int stop=NUMBER_OF_PEPTIDE_ENTRIES_AT_ONCE;
			while (stop<dataAndSourceList.size()) {
				internalWritePeptideQuantLibraryEntriesToConnection(c, inferrer, dataAndSourceList.subList(start, stop), aaConstants);
				if (localizationData.isPresent()) {
					internalWritePeptideLocalizationsToConnection(c, inferrer, localizationData.get(), localizationFDRValues, dataAndSourceList.subList(start, stop));
				}
				start=stop;
				stop=stop+NUMBER_OF_PEPTIDE_ENTRIES_AT_ONCE;
			}
			if (start<dataAndSourceList.size()) {
				internalWritePeptideQuantLibraryEntriesToConnection(c, inferrer, dataAndSourceList.subList(start, dataAndSourceList.size()), aaConstants);
				if (localizationData.isPresent()) {
					internalWritePeptideLocalizationsToConnection(c, inferrer, localizationData.get(), localizationFDRValues, dataAndSourceList.subList(start, dataAndSourceList.size()));
				}
			}

			// FIXME THINK ABOUT WHETHER WE REALLY NEED FRAGMENT LEVEL DATA!
			// PERHAPS WE CAN GET AWAY WITHOUT IT
//			start=0;
//			stop=NUMBER_OF_FRAGMENT_ENTRIES_AT_ONCE;
//			while (stop<dataAndSouceList.size()) {
//				internalWriteFragmentQuantLibraryEntriesToConnection(c, inferrer, dataAndSouceList.subList(start, stop));
//				start=stop;
//				stop=stop+NUMBER_OF_FRAGMENT_ENTRIES_AT_ONCE;
//			}
//			if (start<dataAndSouceList.size()) {
//				internalWriteFragmentQuantLibraryEntriesToConnection(c, inferrer, dataAndSouceList.subList(start, dataAndSouceList.size()));
//			}

			c.commit();
		} finally {
			c.close();
		}
	}

	private static final int NUMBER_OF_PEPTIDE_ENTRIES_AT_ONCE=20;
	// private static final int NUMBER_OF_FRAGMENT_ENTRIES_AT_ONCE=4;

	private void internalWritePeptideQuantLibraryEntriesToConnection(Connection c, Optional<PeakLocationInferrerInterface> inferrer, List<Pair<TransitionRefinementData, String>> dataAndSouceList, AminoAcidConstants aaConstants)
			throws SQLException, IOException {
		int numValidEntries=0;
		for (int i=0; i<dataAndSouceList.size(); i++) {
			if (inferrer.isPresent()) {
				Optional<QuantitativeDIAData> topNIntensity=inferrer.get().getQuantitativeData(dataAndSouceList.get(i).x);
				if (topNIntensity.isPresent()) {
					numValidEntries++;
				}
			} else {
				numValidEntries++;
			}
		}

		StringBuilder peptidePrepString=new StringBuilder("INSERT INTO peptidequants (PrecursorCharge, PeptideModSeq, PeptideSeq, SourceFile, RTInSecondsCenter, "
				+"RTInSecondsStart, RTInSecondsStop, TotalIntensity, NumberOfQuantIons, QuantIonMassLength, "
				+"QuantIonMassArray, QuantIonIntensityLength, QuantIonIntensityArray, BestFragmentCorrelation, BestFragmentDeltaMassPPM, "
				+"MedianChromatogramEncodedLength, MedianChromatogramArray, MedianChromatogramRTEncodedLength, MedianChromatogramRTArray, IdentifiedTICRatio)");
		peptidePrepString.append(" VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
		for (int i=1; i<numValidEntries; i++) {
			peptidePrepString.append(", (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
		}

		PreparedStatement peptidePrep=c.prepareStatement(peptidePrepString.toString());

		try {
			int pepIndex=1;

			for (int i=0; i<dataAndSouceList.size(); i++) {
				Pair<TransitionRefinementData, String> pair=dataAndSouceList.get(i);
				if (inferrer.isPresent()) {
					Optional<QuantitativeDIAData> topNIntensity=inferrer.get().getQuantitativeData(dataAndSouceList.get(i).x);
					if (topNIntensity.isPresent()) {
						pepIndex=prepareQuantData(pair.x, pair.y, inferrer, peptidePrep, pepIndex, aaConstants);
					}
				} else {
					pepIndex=prepareQuantData(pair.x, pair.y, inferrer, peptidePrep, pepIndex, aaConstants);
				}
			}
			if (pepIndex>1) {
				peptidePrep.execute();
			}
		} finally {
			peptidePrep.close();
		}
	}

	private void internalWritePeptideLocalizationsToConnection(Connection c, Optional<PeakLocationInferrerInterface> inferrer, HashMap<String, ModificationLocalizationData> localizationData,
			TObjectDoubleHashMap<String> localizationFDRValues, List<Pair<TransitionRefinementData, String>> dataAndSouceList) throws SQLException, IOException {
		int numStatements=0;
		for (Pair<TransitionRefinementData, String> pair : dataAndSouceList) {
			if (localizationData.containsKey(pair.x.getPeptideModSeq())) {
				numStatements++;
			}
		}

        StringBuilder peptidePrepString=new StringBuilder("INSERT INTO peptidelocalizations (PrecursorCharge, PeptideModSeq, PeptideSeq, SourceFile, LocalizationPeptideModSeq, LocalizationScore, LocalizationFDR, LocalizationIons, NumberOfMods, NumberOfModifiableResidues, IsSiteSpecific, IsLocalized, RTInSecondsCenter, LocalizedIntensity, TotalIntensity)");
        peptidePrepString.append(" VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
        for (int i=1; i<numStatements; i++) {
            peptidePrepString.append(", (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
        }

		PreparedStatement peptidePrep=c.prepareStatement(peptidePrepString.toString());

		try {
			int pepIndex=1;

			for (Pair<TransitionRefinementData, String> pair : dataAndSouceList) {
				pepIndex=prepareLocalizationData(pair.x, pair.y, inferrer, localizationData, localizationFDRValues, peptidePrep, pepIndex);
			}

			if (pepIndex>1) {
				peptidePrep.execute();
			}
		} finally {
			peptidePrep.close();
		}
	}

	public int prepareLocalizationData(TransitionRefinementData data, String sourceFile, Optional<PeakLocationInferrerInterface> inferrer, HashMap<String, ModificationLocalizationData> localizationData,
			TObjectDoubleHashMap<String> localizationFDRValues, PreparedStatement peptidePrep, int index) throws SQLException, IOException {
		ModificationLocalizationData modData=localizationData.get(data.getPeptideModSeq());

		if (modData!=null) {
			peptidePrep.setInt(index++, data.getPrecursorCharge());
			peptidePrep.setString(index++, data.getPeptideModSeq());
			peptidePrep.setString(index++, data.getPeptideSeq());
			peptidePrep.setString(index++, sourceFile);

			peptidePrep.setString(index++, modData.getLocalizationPeptideModSeq().getPeptideAnnotation());
			peptidePrep.setFloat(index++, modData.getLocalizationScore());
			peptidePrep.setDouble(index++, localizationFDRValues.get(data.getPeptideModSeq()));
			peptidePrep.setString(index++, FragmentIon.toArchiveString(modData.getLocalizingIons()));
			peptidePrep.setInt(index++, modData.getNumberOfMods());
			peptidePrep.setInt(index++, modData.getLocalizationPeptideModSeq().getNumModifiableSites());
			peptidePrep.setBoolean(index++, modData.isSiteSpecific());
			peptidePrep.setBoolean(index++, modData.isLocalized());
			peptidePrep.setFloat(index++, modData.getRetentionTimeApexInSeconds());
			peptidePrep.setFloat(index++, modData.getLocalizingIntensity());
			peptidePrep.setFloat(index++, modData.getTotalIntensity());
		}
		return index;
	}

	@SuppressWarnings("unused")
	private void internalWriteFragmentQuantLibraryEntriesToConnection(Connection c, Optional<PeakLocationInferrerInterface> inferrer, List<Pair<TransitionRefinementData, String>> dataAndSouceList)
			throws SQLException, IOException {
		StringBuilder fragmentPrepString=new StringBuilder(
				"INSERT INTO fragmentquants (PrecursorCharge, PeptideModSeq, PeptideSeq, SourceFile, IonType, IonIndex, FragmentMass, Correlation, Background, DeltaMassPPM, Intensity)");
		fragmentPrepString.append(" VALUES (?,?,?,?,?,?,?,?,?,?,?)");
		boolean first=true;
		for (Pair<TransitionRefinementData, String> pair : dataAndSouceList) {
			float[] correlationArray=pair.x.getCorrelationArray();
			for (int i=0; i<correlationArray.length; i++) {
				if (correlationArray[i]>=TransitionRefiner.identificationCorrelationThreshold) {
					if (first) {
						first=false;
					} else {
						fragmentPrepString.append(", (?,?,?,?,?,?,?,?,?,?,?)");
					}
				}
			}
		}

		PreparedStatement fragmentPrep=c.prepareStatement(fragmentPrepString.toString());

		try {
			int fragIndex=1;

			for (Pair<TransitionRefinementData, String> pair : dataAndSouceList) {
				fragIndex=prepareFragmentQuantData(pair.x, pair.y, inferrer, fragmentPrep, fragIndex);
			}
			fragmentPrep.execute();
		} finally {
			fragmentPrep.close();
		}
	}

	public int prepareQuantData(TransitionRefinementData data, String sourceFile, Optional<PeakLocationInferrerInterface> inferrer, PreparedStatement peptidePrep, int index, AminoAcidConstants aaConstants) throws SQLException, IOException {

		QuantitativeDIAData topN;
		if (inferrer.isPresent()) {
			Optional<QuantitativeDIAData> topNIntensity=inferrer.get().getQuantitativeData(data);
			if (!topNIntensity.isPresent()) {
				// not enough ions to quantify peptide
				return index;
			}
			topN=topNIntensity.get();
		} else {
			ArrayList<Peak> peaks=data.getTopNPeaks(TransitionRefiner.quantitativeCorrelationThreshold, Integer.MAX_VALUE);
			Pair<double[], float[]> pair=Peak.toArrays(peaks);
			double[] topNMasses=pair.x;
			float[] topNIntensities=pair.y;
			topN=new QuantitativeDIAData(data.getPeptideModSeq(), data.getPrecursorCharge(), data.getApexRT(), data.getRange(), topNMasses, topNIntensities, aaConstants);
		}

		float[] correlationArray=data.getCorrelationArray();

		Ion[] fragmentMassArray=data.getFragmentMassArray();
		float[] deltaMassArray=data.getDeltaMassArray().get();
		float[] ppmArray=new float[deltaMassArray.length];

		float bestCorrelation=-1.0f;
		float bestDeltaMass=10.0f;
		for (int i=0; i<deltaMassArray.length; i++) {
			ppmArray[i]=deltaMassArray[i]*1000000.0f/(float) fragmentMassArray[i].getMass();
			if (correlationArray[i]>bestCorrelation) {
				bestCorrelation=correlationArray[i];
				bestDeltaMass=ppmArray[i];
			}
		}

		peptidePrep.setInt(index++, topN.getPrecursorCharge());
		peptidePrep.setString(index++, topN.getPeptideModSeq());
		peptidePrep.setString(index++, topN.getPeptideSeq());
		peptidePrep.setString(index++, sourceFile);
		peptidePrep.setFloat(index++, topN.getApexRT());

		peptidePrep.setFloat(index++, topN.getRtScanRange().getStart());
		peptidePrep.setFloat(index++, topN.getRtScanRange().getStop());
		peptidePrep.setFloat(index++, topN.getTIC());
		peptidePrep.setInt(index++, topN.getNumNonZeroPeaks());

		byte[] topNMassesByteArray=ByteConverter.toByteArray(topN.getMassArray());
		peptidePrep.setInt(index++, topNMassesByteArray.length);
		peptidePrep.setBytes(index++, CompressionUtils.compress(topNMassesByteArray));

		byte[] topNIntensitiesByteArray=ByteConverter.toByteArray(topN.getIntensityArray());
		peptidePrep.setInt(index++, topNIntensitiesByteArray.length);
		peptidePrep.setBytes(index++, CompressionUtils.compress(topNIntensitiesByteArray));

		peptidePrep.setFloat(index++, bestCorrelation);
		peptidePrep.setFloat(index++, bestDeltaMass);

		byte[] intensityByteArray=ByteConverter.toByteArray(data.getMedianChromatogram());
		peptidePrep.setInt(index++, intensityByteArray.length);
		peptidePrep.setBytes(index++, CompressionUtils.compress(intensityByteArray));

		if (data.getIdentifiedTICRatio().isPresent()) {
			byte[] rtByteArray=ByteConverter.toByteArray(data.getRtArray().get());
			peptidePrep.setInt(index++, rtByteArray.length);
			peptidePrep.setBytes(index++, CompressionUtils.compress(rtByteArray));

			Float x=data.getIdentifiedTICRatio().get();
			if (x!=null) {
				peptidePrep.setFloat(index++, x);
			} else {
				peptidePrep.setFloat(index++, 0.0f);
			}
		} else {
			peptidePrep.setInt(index++, 0);
			peptidePrep.setBytes(index++, null);
			peptidePrep.setFloat(index++, 0.0f);
		}
		return index;
	}

	public int prepareFragmentQuantData(TransitionRefinementData data, String sourceFile, Optional<PeakLocationInferrerInterface> inferrer, PreparedStatement fragmentPrep, int index)
			throws SQLException, IOException {
		float[] correlationArray=data.getCorrelationArray();
		float[] integrationArray=data.getIntegrationArray();
		float[] backgroundArray=data.getBackgroundArray();

		Ion[] fragmentMassArray=data.getFragmentMassArray();
		float[] deltaMassArray=data.getDeltaMassArray().get();
		float[] ppmArray=new float[deltaMassArray.length];

		float bestCorrelation=-1.0f;
		for (int i=0; i<deltaMassArray.length; i++) {
			ppmArray[i]=deltaMassArray[i]*1000000.0f/(float) fragmentMassArray[i].getMass();
			if (correlationArray[i]>bestCorrelation) {
				bestCorrelation=correlationArray[i];
			}
		}

		for (int i=0; i<correlationArray.length; i++) {
			if (correlationArray[i]>=TransitionRefiner.identificationCorrelationThreshold) {
				fragmentPrep.setInt(index++, data.getPrecursorCharge());
				fragmentPrep.setString(index++, data.getPeptideModSeq());
				fragmentPrep.setString(index++, data.getPeptideSeq());
				fragmentPrep.setString(index++, sourceFile);
				fragmentPrep.setString(index++, IonType.toString(fragmentMassArray[i].getType()));
				fragmentPrep.setInt(index++, fragmentMassArray[i].getIndex());
				fragmentPrep.setDouble(index++, fragmentMassArray[i].getMass());
				fragmentPrep.setFloat(index++, correlationArray[i]);
				fragmentPrep.setFloat(index++, backgroundArray[i]);
				fragmentPrep.setFloat(index++, ppmArray[i]);
				fragmentPrep.setFloat(index++, integrationArray[i]);
			}
		}
		return index;
	}

	public void addTargetDecoyPeptides(ArrayList<PercolatorPeptide> targets, ArrayList<PercolatorPeptide> decoys) throws IOException, SQLException {
		Connection c=getConnection();
		try {
			PreparedStatement prep=c
					.prepareStatement("INSERT INTO peptidescores (PrecursorCharge, PeptideModSeq, PeptideSeq, SourceFile, QValue, PosteriorErrorProbability, IsDecoy) VALUES (?,?,?,?,?,?,?)");
			try {
				for (PercolatorPeptide peptide : targets) {
					prep.setInt(1, peptide.getPrecursorCharge());
					prep.setString(2, peptide.getPeptideModSeq());
					prep.setString(3, peptide.getPeptideSeq());
					prep.setString(4, peptide.getFile());
					prep.setFloat(5, peptide.getQValue());
					prep.setFloat(6, peptide.getPosteriorErrorProb());
					prep.setBoolean(7, peptide.isPSMIDDecoy());
					prep.addBatch();
				}
				for (PercolatorPeptide peptide : decoys) {
					prep.setInt(1, peptide.getPrecursorCharge());
					prep.setString(2, peptide.getPeptideModSeq());
					prep.setString(3, peptide.getPeptideSeq());
					prep.setString(4, peptide.getFile());
					prep.setFloat(5, peptide.getQValue());
					prep.setFloat(6, peptide.getPosteriorErrorProb());
					prep.setBoolean(7, peptide.isPSMIDDecoy());
					prep.addBatch();
				}
				prep.executeBatch();

				c.commit();
			} finally {
				prep.close();
			}
		} finally {
			c.close();
		}
	}

	public ArrayList<ProteinGroupInterface> getProteinGroups() throws IOException, SQLException {
		Connection c=getConnection();
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("SELECT p.proteingroup, p.proteinaccession, p2p.peptideseq, p.qvalue, p.MinimumPeptidePEP FROM proteinscores p, peptidetoprotein p2p WHERE p.proteinaccession=p2p.proteinaccession ORDER BY p.proteingroup");

				ArrayList<ProteinGroupInterface> proteinGroups=new ArrayList<>();

				int previousProteinGroup=-1;
				HashSet<String> accessions=new HashSet<>();
				HashSet<String> peptides=new HashSet<>();
				float prevQValue=0;
				float prevPEP=0;
				while (rs.next()) {
					int proteinGroup=rs.getInt(1);
					String accession=rs.getString(2);
					String peptideseq=rs.getString(3);
					float qvalue=rs.getFloat(4);
					float pep=rs.getFloat(5);

					if (proteinGroup!=previousProteinGroup) {
						previousProteinGroup=proteinGroup;
						if (accessions.size()>0) {
							proteinGroups.add(new PercolatorProteinGroup(accessions.toArray(new String[accessions.size()]), peptides.toArray(new String[peptides.size()]), prevQValue, prevPEP));
						}
						accessions.clear();
						accessions.add(accession);
						peptides.clear();
						peptides.add(peptideseq);
						prevQValue=qvalue;
						prevPEP=pep;
					} else {
						accessions.add(accession);
						peptides.add(peptideseq);
					}
				}
				if (accessions.size()>0) {
					proteinGroups.add(new PercolatorProteinGroup(accessions.toArray(new String[accessions.size()]), peptides.toArray(new String[peptides.size()]), prevQValue, prevPEP));
				}

				return proteinGroups;
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	public void addTargetDecoyProteins(String fileID, ArrayList<PercolatorProteinGroup> targets, ArrayList<PercolatorProteinGroup> decoys) throws IOException, SQLException {
		Connection c=getConnection();
		try {
			PreparedStatement prep=c
					.prepareStatement("INSERT INTO proteinscores (ProteinGroup, ProteinAccession, SourceFile, QValue, MinimumPeptidePEP, IsDecoy) VALUES (?,?,?,?,?,?)");
			try {
				int group=0;
				for (PercolatorProteinGroup protein : targets) {
					group++;
					for (String accession : protein.getAccessions()) {
						prep.setInt(1, group);
						prep.setString(2, accession);
						prep.setString(3, fileID);
						prep.setFloat(4, protein.getQValue());
						prep.setFloat(5, protein.getPosteriorErrorProb());
						prep.setBoolean(6, false);
						prep.addBatch();
					}
				}
				for (PercolatorProteinGroup protein : decoys) {
					group++;
					for (String accession : protein.getAccessions()) {
						prep.setInt(1, group);
						prep.setString(2, accession);
						prep.setString(3, fileID);
						prep.setFloat(4, protein.getQValue());
						prep.setFloat(5, protein.getPosteriorErrorProb());
						prep.setBoolean(6, true);
						prep.addBatch();
					}
				}
				prep.executeBatch();

				c.commit();
			} finally {
				prep.close();
			}
		} finally {
			c.close();
		}
	}

	public void addProteinsFromPercolator(ArrayList<PercolatorPeptide> entries) throws IOException, SQLException {
		Connection c=getConnection();
		try {
			PreparedStatement proteinPrep=c.prepareStatement("INSERT OR IGNORE INTO peptidetoprotein (PeptideSeq, isDecoy, ProteinAccession) VALUES (?,?,?)");
			try {
				for (PercolatorPeptide percolatorPeptide : entries) {
					proteinPrep.setString(1, percolatorPeptide.getPeptideSeq());
					proteinPrep.setBoolean(2, percolatorPeptide.isPSMIDDecoy());
					for (String acc : PSMData.stringToAccessions(percolatorPeptide.getProteinIDs())) {
						proteinPrep.setString(3, acc);
						proteinPrep.addBatch();
					}
				}
				proteinPrep.executeBatch();

				c.commit();
			} finally {
				proteinPrep.close();
			}
		} finally {
			c.close();
		}
	}

	public void addProteinsFromEntries(ArrayList<LibraryEntry> entries) throws IOException, SQLException {
		Connection c=getConnection();
		try {
			PreparedStatement proteinPrep=c.prepareStatement("INSERT OR IGNORE INTO peptidetoprotein (PeptideSeq, IsDecoy, ProteinAccession) VALUES (?,?,?)");
			try {
				HashMap<String, HashSet<String>> targetAccessionsByPeptide=new HashMap<>();
				HashMap<String, HashSet<String>> decoyAccessionsByPeptide=new HashMap<>();

				for (LibraryEntry entry : entries) {
					HashMap<String, HashSet<String>> map;
					if (entry.isDecoy()) {
						map=decoyAccessionsByPeptide;
					} else {
						map=targetAccessionsByPeptide;
					}
					
					HashSet<String> accessions=map.get(entry.getPeptideSeq());
					if (accessions==null) {
						accessions=new HashSet<>();
						map.put(entry.getPeptideSeq(), accessions);
					}
					accessions.addAll(entry.getAccessions());
				}
				
				for (Entry<String, HashSet<String>> entry : targetAccessionsByPeptide.entrySet()) {
					proteinPrep.setString(1, entry.getKey());
					proteinPrep.setBoolean(2, false);
					for (String acc : entry.getValue()) {
						proteinPrep.setString(3, acc);
						proteinPrep.addBatch();
					}
				}
				
				for (Entry<String, HashSet<String>> entry : decoyAccessionsByPeptide.entrySet()) {
					proteinPrep.setString(1, entry.getKey());
					proteinPrep.setBoolean(2, true);
					for (String acc : entry.getValue()) {
						proteinPrep.setString(3, acc);
						proteinPrep.addBatch();
					}
				}
				
				proteinPrep.executeBatch();

				c.commit();
			} finally {
				proteinPrep.close();
			}
		} finally {
			c.close();
		}
	}

	public void addEntries(ArrayList<LibraryEntry> entries) throws IOException, SQLException {
		Connection c=getConnection();
		try {
			PreparedStatement prep=c.prepareStatement(
					"INSERT INTO entries (PrecursorMZ, PrecursorCharge, PeptideModSeq, PeptideSeq, Copies, RTInSeconds, Score, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray, CorrelationEncodedLength, CorrelationArray, RTInSecondsStart, RTInSecondsStop, MedianChromatogramEncodedLength, MedianChromatogramArray, SourceFile) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
			try {
				for (LibraryEntry entry : entries) {
					if (entry.getAccessions().size()==0)
						continue;
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

					if (entry.getMassArray().length!=entry.getIntensityArray().length) {
						throw new EncyclopediaException("Mass/Intensity array length mismatch! "+entry.getMassArray().length+" != "+entry.getIntensityArray().length+" FOR "+entry.getPeptideModSeq());
					}

					if (entry instanceof Chromatogram) {
						Chromatogram cast=(Chromatogram) entry;

						byte[] correlationByteArray=ByteConverter.toByteArray(cast.getCorrelationArray());
						if (entry.getMassArray().length!=entry.getIntensityArray().length) {
							throw new EncyclopediaException(
									"Mass/Correlation array length mismatch! "+entry.getMassArray().length+" != "+entry.getIntensityArray().length+" FOR "+entry.getPeptideModSeq());
						}

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
				}
				prep.executeBatch();

				c.commit();
			} finally {
				prep.close();
			}
		} finally {
			c.close();
		}
	}

	public Connection getConnection() throws IOException {
		return getConnection(tempFile);
	}

//	public HashMap<PeptidePrecursor, ArrayList<LibraryEntry>> getEntries(ArrayList<PeptidePrecursor> entries, boolean sqrt) throws IOException, SQLException, DataFormatException {
//		HashMap<PeptidePrecursor, ArrayList<LibraryEntry>> map=new HashMap<PeptidePrecursor, ArrayList<LibraryEntry>>();
//
//		try (Connection c=getConnection()) {
//			try (PreparedStatement prep=c.prepareStatement("select "+"e.PrecursorMZ, "+"e.PrecursorCharge, "+"e.PeptideModSeq, "+"e.Copies, "+"e.RTInSeconds, "+"e.Score, "+"e.MassEncodedLength, "
//					+"e.MassArray, "+"e.IntensityEncodedLength, "+"e.IntensityArray, "+"e.CorrelationEncodedLength, "+"e.CorrelationArray blob, "+"e.RTInSecondsStart, "+"e.RTInSecondsStop,"
//					+"e.MedianChromatogramEncodedLength, "+"e.MedianChromatogramArray, "+"group_concat(p.ProteinAccession, '"+PSMData.ACCESSION_TOKEN+"') ProteinAccessions, "+"e.SourceFile "+"from "
//					+"entries e "+"left join peptidetoprotein p "+"on "+"e.PeptideSeq=p.PeptideSeq "+"and not p.isdecoy "+"where e.PeptideModSeq = ? "+"and e.PrecursorCharge = ? "+"group by e.rowid;")) {
//				for (PeptidePrecursor precursor : entries) {
//					prep.setString(1, precursor.getPeptideModSeq());
//					prep.setByte(2, precursor.getPrecursorCharge());
//
//					ResultSet rs=prep.executeQuery();
//					ArrayList<LibraryEntry> entry=extractEntries(sqrt, rs);
//					map.put(precursor, entry);
//				}
//			}
//		}
//
//		return map;
//	}

	/*
	 * (non-Javadoc)
	 *
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface#
	 * getEntries(edu.washington.gs.maccoss.encyclopedia.datastructures.Range)
	 */
	@Override
	public ArrayList<LibraryEntry> getEntries(String peptideModSeq, byte charge, boolean sqrt) throws IOException, SQLException, DataFormatException {
		try (Connection c=getConnection()) {
			try (PreparedStatement s=c.prepareStatement("select "+"e.PrecursorMZ, "+"e.PrecursorCharge, "+"e.PeptideModSeq, "+"e.Copies, "+"e.RTInSeconds, "+"e.Score, "+"e.MassEncodedLength, "
					+"e.MassArray, "+"e.IntensityEncodedLength, "+"e.IntensityArray, "+"e.CorrelationEncodedLength, "+"e.CorrelationArray blob, "+"e.RTInSecondsStart, "+"e.RTInSecondsStop, "
					+"e.MedianChromatogramEncodedLength, "+"e.MedianChromatogramArray, "+"group_concat(p.ProteinAccession, '"+PSMData.ACCESSION_TOKEN+"') ProteinAccessions, "+"e.SourceFile "+"from "
					+"entries e "+"left join peptidetoprotein p "+"on "+"e.PeptideSeq=p.PeptideSeq "+"and not p.isdecoy "+"where e.PeptideModSeq = ? "+"and e.PrecursorCharge = ? "+"group by e.rowid;")) {
				s.setString(1, peptideModSeq);
				s.setByte(2, charge);
				ResultSet rs=s.executeQuery();

				// Don't bother with modification mass munging
				final ArrayList<LibraryEntry> libraryEntries = extractEntries(sqrt, rs, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()));

				//TODO: what if there's an entry that would have its seq munged to match the query!?!?

				// Check that mod masses weren't affected (the results have the mod seq that was requested)
				if (!libraryEntries.stream().map(LibraryEntry::getPeptideModSeq).allMatch(Predicate.isEqual(peptideModSeq))) {
					throw new IllegalStateException("Result had mismatched modSeq!");
				}

				return libraryEntries;
			}
		}
	}

	private ArrayList<LibraryEntry> extractEntries(boolean sqrt, ResultSet rs, AminoAcidConstants aaConstants) throws SQLException, IOException, DataFormatException {
		String peptideModSeq;
		ArrayList<LibraryEntry> entry=new ArrayList<LibraryEntry>();
		while (rs.next()) {

			double precursorMZ=rs.getDouble(1);
			byte precursorCharge=(byte) rs.getInt(2);
			peptideModSeq=PeptideUtils.getCorrectedMasses(rs.getString(3), aaConstants);
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

			if (massArray.length!=intensityArray.length) {
				throw new EncyclopediaException("Mass/Intensity array length mismatch! "+massArray.length+" != "+intensityArray.length+" FOR "+peptideModSeq);
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
				if (massArray.length!=correlationArray.length) {
					throw new EncyclopediaException("Mass/Correlation array length mismatch! "+massArray.length+" != "+correlationArray.length+" FOR "+peptideModSeq);
				}

				rtInSecondsStart=rs.getFloat(13);
				rtInSecondsStop=rs.getFloat(14);
				int medianChromatogramEncodedLength=rs.getInt(15);
				medianChromatogramArray=ByteConverter.toFloatArray(CompressionUtils.decompress(rs.getBytes(16), medianChromatogramEncodedLength));
			}

			HashSet<String> accessions=PSMData.stringToAccessions(rs.getString(17));
			String sourceFile=rs.getString(18);
			if (correlationEncodedLength==0) {
				entry.add(new LibraryEntry(sourceFile, accessions, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, aaConstants));
			} else {
				entry.add(new ChromatogramLibraryEntry(sourceFile, accessions, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray,
						correlationArray, medianChromatogramArray, new Range(rtInSecondsStart, rtInSecondsStop), aaConstants));
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
					return new Range((float) min, (float) max);
				}

				return new Range(0, 0);
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}
	
	public HashMap<String, String> getAccessions(Collection<String> peptideSeqs) throws IOException, SQLException, DataFormatException {
		HashMap<String, String> accessions=new HashMap<>();
		
		try (Connection c = getConnection()) {	
			try (PreparedStatement s = c.prepareStatement("select group_concat(ProteinAccession, '"+PSMData.ACCESSION_TOKEN+"') ProteinAccessions from peptidetoprotein where peptideseq = ?")) {
				for (String peptideseq : peptideSeqs) {
					s.setString(1, peptideseq);
					ResultSet rs=s.executeQuery();
					while (rs.next()) {
						accessions.put(peptideseq, rs.getString(1));
					}
				}

			} finally {
				c.close();
			}
		}
		return accessions;
	}
	
	PreparedStatement accessionsStatement=null;
	
	public HashSet<String> getAccessions(String peptideSeq) throws IOException, SQLException, DataFormatException {
		Connection connection = getConnection();
		Statement createStatement = connection.createStatement();
		ResultSet rs=createStatement.executeQuery("select proteinaccession from peptidetoprotein where peptideseq = '"+peptideSeq+"'");
		HashSet<String> accessions=new HashSet<>();
		while (rs.next()) {
			accessions.add(rs.getString(1));
		}
		rs.close();
		createStatement.close();
		connection.close();
		return accessions;
	}
	
	@Override
	public ArrayList<LibraryEntry> getUnlinkedEntries(Range precursorMz, boolean sqrt, AminoAcidConstants aaConstants) throws IOException, SQLException, DataFormatException {
		try (Connection c = getConnection()) {
			try (PreparedStatement s = c.prepareStatement("select " 
					+ "e.PrecursorMZ, " 
					+ "e.PrecursorCharge, "
					+ "e.PeptideModSeq, " 
					+ "e.PeptideSeq, "
					+ "e.Copies, " 
					+ "e.RTInSeconds, " 
					+ "e.Score, " 
					+ "e.MassEncodedLength, "
					+ "e.MassArray, " 
					+ "e.IntensityEncodedLength, " 
					+ "e.IntensityArray, "
					+ "e.SourceFile " 
					+ "from " 
					+ "entries e "
					+ "where e.PrecursorMz between ? and ?;")) {
				s.setFloat(1, precursorMz.getStart());
				s.setFloat(2, precursorMz.getStop());
				ResultSet rs=s.executeQuery();

				ArrayList<LibraryEntry> entry=new ArrayList<LibraryEntry>();
				while (rs.next()) {

					double precursorMZ=rs.getDouble(1);
					byte precursorCharge=(byte) rs.getInt(2);
					String peptideModSeq=PeptideUtils.getCorrectedMasses(rs.getString(3), aaConstants);
					String peptideSeq=rs.getString(4);
					int copies=rs.getInt(5);
					float retentionTime=rs.getFloat(6);
					float score=rs.getFloat(7);
					int massEncodedLength=rs.getInt(8);
					double[] massArray=ByteConverter.toDoubleArray(CompressionUtils.decompress(rs.getBytes(9), massEncodedLength));
					int intensityEncodedLength=rs.getInt(10);
					float[] intensityArray=ByteConverter.toFloatArray(CompressionUtils.decompress(rs.getBytes(11), intensityEncodedLength));
					if (sqrt) {
						intensityArray=General.protectedSqrt(intensityArray);
					}
					String sourceFile=rs.getString(12);
					entry.add(new UnlinkedLibraryEntry(sourceFile, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, aaConstants, false, false, false, peptideSeq, this));
				}

				return entry;
			} finally {
				c.close();
			}
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface#
	 * getEntries(edu.washington.gs.maccoss.encyclopedia.datastructures.Range)
	 */
	@Override
	public ArrayList<LibraryEntry> getEntries(Range precursorMz, boolean sqrt, AminoAcidConstants aaConstants) throws IOException, SQLException, DataFormatException {
		try (Connection c=getConnection()) {
			try (PreparedStatement s=c.prepareStatement("select "+"e.PrecursorMZ, "+"e.PrecursorCharge, "+"e.PeptideModSeq, "+"e.Copies, "+"e.RTInSeconds, "+"e.Score, "+"e.MassEncodedLength, "
					+"e.MassArray, "+"e.IntensityEncodedLength, "+"e.IntensityArray, "+"e.CorrelationEncodedLength, "+"e.CorrelationArray blob, "+"e.RTInSecondsStart, "+"e.RTInSecondsStop, "
					+"e.MedianChromatogramEncodedLength, "+"e.MedianChromatogramArray, "+"group_concat(p.ProteinAccession, '"+PSMData.ACCESSION_TOKEN+"') ProteinAccessions, "+"e.SourceFile "+"from "
					+"entries e "+"left join peptidetoprotein p "+"on "+"e.PeptideSeq=p.PeptideSeq "+"and not p.isdecoy "+"where e.PrecursorMz between ? and ? "+"group by e.rowid;")) {
				s.setFloat(1, precursorMz.getStart());
				s.setFloat(2, precursorMz.getStop());
				ResultSet rs=s.executeQuery();

				ArrayList<LibraryEntry> entry=new ArrayList<LibraryEntry>();
				while (rs.next()) {

					double precursorMZ=rs.getDouble(1);
					byte precursorCharge=(byte) rs.getInt(2);
					String peptideModSeq=PeptideUtils.getCorrectedMasses(rs.getString(3), aaConstants);
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
						entry.add(new LibraryEntry(sourceFile, accessions, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, aaConstants));
					} else {
						entry.add(new ChromatogramLibraryEntry(sourceFile, accessions, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray,
								correlationArray, medianChromatogramArray, new Range(rtInSecondsStart, rtInSecondsStop), aaConstants));
					}
				}

				return entry;
			}
		}
	}

	public ArrayList<PeptidePrecursor> getAllPeptidePrecursors(AminoAcidConstants aaConstants) throws IOException, SQLException, DataFormatException {
		try (Connection c=getConnection()) {
			try (PreparedStatement s=c.prepareStatement("select "+"e.PeptideModSeq, "+"e.PrecursorCharge "+"from "
					+"entries e "+"left join peptidetoprotein p "+"on "+"e.PeptideSeq=p.PeptideSeq "+"and not p.isdecoy "+"group by e.rowid")) {

				ResultSet rs=s.executeQuery();

				ArrayList<PeptidePrecursor> entry=new ArrayList<PeptidePrecursor>();
				while (rs.next()) {
					String peptideModSeq=PeptideUtils.getCorrectedMasses(rs.getString(1), aaConstants);
					byte precursorCharge=(byte) rs.getInt(2);
					SimplePeptidePrecursor peptide=new SimplePeptidePrecursor(peptideModSeq, precursorCharge, aaConstants);
					entry.add(peptide);
				}

				return entry;
			}
		}
	}

	public ArrayList<LibraryEntry> getAllEntries(boolean sqrt, AminoAcidConstants aaConstants) throws IOException, SQLException, DataFormatException {
		try (Connection c=getConnection()) {
			try (PreparedStatement s=c.prepareStatement("select "+"e.PrecursorMZ, "+"e.PrecursorCharge, "+"e.PeptideModSeq, "+"e.Copies, "+"e.RTInSeconds, "+"e.Score, "+"e.MassEncodedLength, "
					+"e.MassArray, "+"e.IntensityEncodedLength, "+"e.IntensityArray, "+"e.CorrelationEncodedLength, "+"e.CorrelationArray blob, "+"e.RTInSecondsStart, "+"e.RTInSecondsStop, "
					+"e.MedianChromatogramEncodedLength, "+"e.MedianChromatogramArray, "+"group_concat(p.ProteinAccession, '"+PSMData.ACCESSION_TOKEN+"') ProteinAccessions, "+"e.SourceFile "+"from "
					+"entries e "+"left join peptidetoprotein p "+"on "+"e.PeptideSeq=p.PeptideSeq "+"and not p.isdecoy "+"group by e.rowid")) {

				ResultSet rs=s.executeQuery();

				ArrayList<LibraryEntry> entry=new ArrayList<LibraryEntry>();
				while (rs.next()) {

					double precursorMZ=rs.getDouble(1);
					byte precursorCharge=(byte) rs.getInt(2);
					String peptideModSeq=PeptideUtils.getCorrectedMasses(rs.getString(3), aaConstants);
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
						entry.add(new LibraryEntry(sourceFile, accessions, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, aaConstants));
					} else {
						entry.add(new ChromatogramLibraryEntry(sourceFile, accessions, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray,
								correlationArray, medianChromatogramArray, new Range(rtInSecondsStart, rtInSecondsStop), aaConstants));
					}
				}

				return entry;
			}
		}
	}

    public ArrayList<LocalizedLibraryEntry> getAllLocalizedEntries(float minimumLocalizationScore, boolean requireLocalization, PeptideModification mod, boolean sqrt, AminoAcidConstants aaConstants) throws IOException, SQLException, DataFormatException {
        try (Connection c = getConnection()) {
            String sql = "select " +
                "e.PrecursorMZ, " +
                "e.PrecursorCharge, " +
                "e.PeptideModSeq, " +
                "e.Copies, " +
                "l.RTInSecondsCenter, " +
                "e.Score, " +
                "e.MassEncodedLength, " +
                "e.MassArray, " +
                "e.IntensityEncodedLength, " +
                "e.IntensityArray, " +
                "e.CorrelationEncodedLength, " +
                "e.CorrelationArray blob, " +
                "e.RTInSecondsStart, " +
                "e.RTInSecondsStop, " +
                "e.MedianChromatogramEncodedLength, " +
                "e.MedianChromatogramArray, " +
                "group_concat(p.ProteinAccession, '" + PSMData.ACCESSION_TOKEN + "') ProteinAccessions, " +
                "e.SourceFile, " +
                "l.LocalizationPeptideModSeq, " +
                "l.LocalizationScore, " +
                "l.LocalizationIons, " +
                "l.NumberOfMods, " +
                "l.NumberOfModifiableResidues, " +
                "l.isSiteSpecific "+
                "from " +
                "peptidelocalizations l, " +
                "entries e " +
                "left join peptidetoprotein p " +
                "on "+"e.PeptideSeq=p.PeptideSeq and not p.isdecoy " +
                "where " +
                (requireLocalization?"l.isLocalized=1 and ":"")+
                "l.LocalizationScore>="+minimumLocalizationScore+" and " +
                "e.PeptideModSeq=l.PeptideModSeq and " +
                "e.PrecursorCharge=l.PrecursorCharge and " +
                "e.SourceFile=l.SourceFile and " +
                "e.PeptideSeq=p.PeptideSeq " +
                "group by e.rowid";
            try (PreparedStatement s = c.prepareStatement(
                sql
            )) {
                ResultSet rs = s.executeQuery();

				ArrayList<LocalizedLibraryEntry> entry=new ArrayList<LocalizedLibraryEntry>();
				while (rs.next()) {

					double precursorMZ=rs.getDouble(1);
					byte precursorCharge=(byte) rs.getInt(2);
					String peptideModSeq=PeptideUtils.getCorrectedMasses(rs.getString(3), aaConstants);

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

					AmbiguousPeptideModSeq peptideAnnotation=AmbiguousPeptideModSeq.getAmbiguousPeptideModSeq(rs.getString(19), mod, aaConstants);
					float localizationScore=rs.getFloat(20);
					FragmentIon[] localizationIons=FragmentIon.fromArchiveString(rs.getString(21));
					int numberOfMods=rs.getInt(22);
					int numberOfModifiableResidues=rs.getInt(23);
					boolean isSiteSpecific=rs.getBoolean(24);

					if (correlationEncodedLength!=0) {
						entry.add(new LocalizedLibraryEntry(sourceFile, accessions, 1, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray,
								correlationArray, medianChromatogramArray, new Range(rtInSecondsStart, rtInSecondsStop), peptideAnnotation, localizationScore, localizationIons,
								numberOfModifiableResidues, numberOfMods, isSiteSpecific, aaConstants));
					}
				}

				return entry;
			}
		}
	}

	private static void populatePeptideToProtein(Connection c) throws SQLException {
		// TODO: if you re-enable this, be sure to close these statements which are not being closed
//		System.out.println("POPULATING PEPTIDE TO PROTEIN");
//		System.out.println("entries SIZE: "+c.prepareStatement("select count(*) from entries").executeQuery().getInt(1));
//		System.out.println("proteins SIZE: "+c.prepareStatement("select count(*) from proteins").executeQuery().getInt(1));
//		System.out.println("peptidetoprotein SIZE: "+c.prepareStatement("select count(*) from peptidetoprotein").executeQuery().getInt(1));

		try (PreparedStatement s=c.prepareStatement("select * from proteins p;")) {
			s.execute();
			try (ResultSet rs=s.getResultSet()) {
				try (PreparedStatement ins=c.prepareStatement("insert into peptidetoprotein (peptideseq, isDecoy, proteinaccession) values (?, ?, ?);")) {
					int i=0;
					while (rs.next()) {
						ins.setString(1, rs.getString("peptideseq"));
						for (String acc : PSMData.stringToAccessions(rs.getString("proteinaccessions"))) {
							ins.setBoolean(2, false); // everything that was
														// originally added to
														// proteins is not a
														// decoy
							ins.setString(3, acc);
							ins.addBatch();
							i++;
						}
						if (i>1024) {
							// only execute batches outside inner loop to avoid
							// losing the peptide sequence param
							ins.executeBatch();
							ins.clearBatch();
							ins.clearParameters();
							i=0;
						}
					}
					ins.executeBatch();
				} catch (Exception e) {
					e.printStackTrace();
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		// TODO: if you re-enable this, be sure to close these statements which are not being closed
//		System.out.println("CHECKING PEPTIDE TO PROTEIN");
//		System.out.println("entries SIZE: "+c.prepareStatement("select count(*) from entries").executeQuery().getInt(1));
//		System.out.println("proteins SIZE: "+c.prepareStatement("select count(*) from proteins").executeQuery().getInt(1));
//		System.out.println("peptidetoprotein SIZE: "+c.prepareStatement("select count(*) from peptidetoprotein").executeQuery().getInt(1));

	}

	private void createNewTables() throws IOException, SQLException {
		Connection c=getConnection();
		try {
			Statement s=c.createStatement();
			try {
				// lack of a metadata table implies this is a new file. no patches are needed, and all tables will
				// be created
				Version version = doesTableExist(c, "metadata") ? getVersion() : null;

				if (version!=null) {
					if (userFile!=null) {
						Logger.logLine("Opening library "+userFile.getName()+" (version: "+version+")");
					}

					if (new Version(0, 1, 2).amIAbove(version) && version.amIAbove(new Version(0, 0, 9))) {
						if (userFile!=null) {
							Logger.logLine("Updating library to "+new Version(0, 1, 2));
						}
						s.execute("ALTER TABLE entries ADD COLUMN CorrelationEncodedLength int");
						s.execute("ALTER TABLE entries ADD COLUMN CorrelationArray blob");
						s.execute("ALTER TABLE entries ADD COLUMN RTInSecondsStart double");
						s.execute("ALTER TABLE entries ADD COLUMN RTInSecondsStop double");
						s.execute("ALTER TABLE entries ADD COLUMN MedianChromatogramEncodedLength int");
						s.execute("ALTER TABLE entries ADD COLUMN MedianChromatogramArray blob");
					}
					if (new Version(0, 1, 4).amIAbove(version)&&version.amIAbove(new Version(0, 1, 2))) {
						if (userFile!=null) {
							Logger.logLine("Updating library to "+new Version(0, 1, 4));
						}
						s.execute("ALTER TABLE fragmentquants ADD COLUMN Background double");
						s.execute("ALTER TABLE fragmentquants ADD COLUMN PeptideSeq string");
						s.execute("ALTER TABLE peptidequants ADD COLUMN PeptideSeq string");
					}

					if (new Version(0, 1, 5).amIAbove(version)&&version.amIAbove(new Version(0, 1, 2))) {
						if (userFile!=null) {
							Logger.logLine("Updating library to "+new Version(0, 1, 5));
						}
						s.execute("ALTER TABLE peptidequants ADD COLUMN RTInSecondsCenter double");
					}

					if (new Version(0, 1, 6).amIAbove(version)&&version.amIAbove(new Version(0, 1, 2))) {
						if (userFile!=null) {
							Logger.logLine("Updating library to "+new Version(0, 1, 6));
						}
						s.execute("ALTER TABLE peptidequants ADD COLUMN IdentifiedTICRatio double");
					}

					if (new Version(0, 1, 7).amIAbove(version)&&version.amIAbove(new Version(0, 1, 2))) {
						if (userFile!=null) {
							Logger.logLine("Updating library to "+new Version(0, 1, 7));
						}
						s.execute("ALTER TABLE fragmentquants ADD COLUMN IonIndex int");
					}

					if (new Version(0, 1, 8).amIAbove(version)&&version.amIAbove(new Version(0, 1, 2))) {
						if (userFile!=null) {
							Logger.logLine("Updating library to "+new Version(0, 1, 8));
						}
						s.execute("ALTER TABLE peptidequants ADD COLUMN QuantIonMassLength int");
						s.execute("ALTER TABLE peptidequants ADD COLUMN QuantIonMassArray blob");
					}

					if (new Version(0, 1, 10).amIAbove(version)&&version.amIAbove(new Version(0, 1, 2))) {
						if (userFile!=null) {
							Logger.logLine("Updating library to "+new Version(0, 1, 10));
						}
						s.execute("CREATE TABLE IF NOT EXISTS peptidetoprotein ("+"PeptideSeq string not null,"+"ProteinAccession string not null"+");");
						s.execute("ALTER TABLE peptidetoprotein ADD COLUMN isDecoy boolean");

						populatePeptideToProtein(c);
						s.execute("DROP TABLE proteins;");
					}

					if (new Version(0, 1, 11).amIAbove(version)&&version.amIAbove(new Version(0, 1, 2))) {
						if (userFile!=null) {
							Logger.logLine("Updating library to "+new Version(0, 1, 11));
						}
						s.execute("ALTER TABLE peptidequants ADD COLUMN QuantIonIntensityLength int");
						s.execute("ALTER TABLE peptidequants ADD COLUMN QuantIonIntensityArray blob");
						s.execute("ALTER TABLE peptidequants ADD COLUMN MedianChromatogramRTEncodedLength int");
						s.execute("ALTER TABLE peptidequants ADD COLUMN MedianChromatogramRTArray blob");
					}

					if (new Version(0, 1, 12).amIAbove(version)&&version.amIAbove(new Version(0, 1, 9))) {
						if (userFile!=null) {
							Logger.logLine("Updating library to "+new Version(0, 1, 12));
						}
						s.execute("ALTER TABLE peptidetoprotein ADD COLUMN isDecoy boolean");
					}
				}

				// UNIQUE constraints cost as much as an index and can't
				// add/drop them, so we have to live without the constraint and
				// deal with it in code
				s.execute("CREATE TABLE IF NOT EXISTS metadata ( "+"Key string not null, Value string not null "+")"); // +"UNIQUE
																														// (Key)
																														// )");

				s.execute("CREATE TABLE IF NOT EXISTS entries ( "
						+"PrecursorMz double not null, PrecursorCharge int not null, PeptideModSeq string not null, PeptideSeq string not null, Copies int not null, RTInSeconds double not null, Score double not null, MassEncodedLength int not null, MassArray blob not null, IntensityEncodedLength int not null, IntensityArray blob not null, CorrelationEncodedLength int, CorrelationArray blob, RTInSecondsStart double, RTInSecondsStop double, MedianChromatogramEncodedLength int, MedianChromatogramArray blob, SourceFile string not null "
						+")"); // +"UNIQUE (PrecursorCharge, PeptideModSeq,
								// SourceFile) )");

				s.execute("CREATE TABLE IF NOT EXISTS peptidetoprotein ("+"PeptideSeq string not null,"+"isDecoy boolean,"+"ProteinAccession string not null"+");");

				s.execute("CREATE TABLE IF NOT EXISTS peptidequants ( "
						+"PrecursorCharge int not null, PeptideModSeq string not null, PeptideSeq string not null, SourceFile string not null, RTInSecondsCenter double not null, RTInSecondsStart double not null, RTInSecondsStop double not null, TotalIntensity double not null, NumberOfQuantIons int not null, QuantIonMassLength int not null, QuantIonMassArray blob not null, QuantIonIntensityLength int, QuantIonIntensityArray blob, BestFragmentCorrelation double not null, BestFragmentDeltaMassPPM double not null, MedianChromatogramEncodedLength int not null, MedianChromatogramArray blob not null, MedianChromatogramRTEncodedLength int, MedianChromatogramRTArray blob, IdentifiedTICRatio double not null "
						+")"); // +"UNIQUE (PrecursorCharge, PeptideModSeq,
								// SourceFile) )");

				s.execute("CREATE TABLE IF NOT EXISTS peptidelocalizations ( "
						+"PrecursorCharge int not null, PeptideModSeq string not null, PeptideSeq string not null, SourceFile string not null, LocalizationPeptideModSeq string, LocalizationScore double, LocalizationFDR double, LocalizationIons string, NumberOfMods int, NumberOfModifiableResidues int, IsSiteSpecific boolean, IsLocalized boolean, RTInSecondsCenter double, LocalizedIntensity double, TotalIntensity double "
						+")"); // +"UNIQUE (PrecursorCharge, PeptideModSeq,
								// SourceFile) )");

				s.execute("CREATE TABLE IF NOT EXISTS fragmentquants ( "
						+"PrecursorCharge int not null, PeptideModSeq string not null, PeptideSeq string not null, SourceFile string not null, IonType string not null, IonIndex int not null, FragmentMass double not null, Correlation double not null, Background double not null, DeltaMassPPM double not null, Intensity double not null "
						+")");

				s.execute("CREATE TABLE IF NOT EXISTS peptidescores ( "
						+"PrecursorCharge int not null, PeptideModSeq string not null, PeptideSeq string not null, SourceFile string not null, QValue double not null, PosteriorErrorProbability double not null, IsDecoy boolean not null "
						+")"); // +"UNIQUE (PrecursorCharge, PeptideModSeq,
								// SourceFile) )");

				s.execute("CREATE TABLE IF NOT EXISTS proteinscores ( "
						+"ProteinGroup int not null, ProteinAccession string not null, SourceFile string not null, QValue double not null, MinimumPeptidePEP double not null, IsDecoy boolean not null "
						+")");

				s.execute("CREATE TABLE IF NOT EXISTS retentiontimes (SourceFile string not null, Library float not null, Actual float not null, Predicted float not null, Delta float not null, Probability float not null, Decoy boolean, PeptideModSeq string)");

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
				s.execute("drop index if exists 'Key_Metadata_index'");

				s.execute("drop index if exists 'PeptideModSeq_PrecursorCharge_SourceFile_Entries_index'");
				s.execute("drop index if exists 'PeptideSeq_Entries_index'");
				s.execute("drop index if exists 'PrecursorMz_Entries_index'");

				s.execute("drop index if exists 'PeptideModSeq_PrecursorCharge_SourceFile_Peptides_index'");
				s.execute("drop index if exists 'PeptideSeq_Peptides_index'");

				s.execute("drop index if exists 'PeptideModSeq_PrecursorCharge_SourceFile_Localizations_index'");
				s.execute("drop index if exists 'PeptideSeq_Localizations_index'");

				s.execute("drop index if exists 'PeptideModSeq_PrecursorCharge_SourceFile_Fragments_index'");
				s.execute("drop index if exists 'PeptideSeq_Fragments_index'");

				s.execute("drop index if exists 'ProteinGroup_ProteinScores_index'");
				s.execute("drop index if exists 'ProteinAccession_ProteinScores_index'");

				s.execute("drop index if exists 'PeptideModSeq_PrecursorCharge_SourceFile_Scores_index'");
				s.execute("drop index if exists 'PeptideSeq_Scores_index'");

				s.execute("drop index if exists 'ProteinAccession_PeptideToProtein_index'");
				s.execute("drop index if exists 'PeptideSeq_PeptideToProtein_index'");

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
				s.execute("create index if not exists 'Key_Metadata_index' on 'metadata' ('Key' ASC)");

				s.execute("create index if not exists 'PeptideModSeq_PrecursorCharge_SourceFile_Entries_index' on 'entries' ('PeptideModSeq' ASC, 'PrecursorCharge' ASC, 'SourceFile' ASC)");
				s.execute("create index if not exists 'PeptideSeq_Entries_index' on 'entries' ('PeptideSeq' ASC)");
				s.execute("create index if not exists 'PrecursorMz_Entries_index' on 'entries' ('PrecursorMz' ASC)");

				s.execute("create index if not exists 'PeptideModSeq_PrecursorCharge_SourceFile_Peptides_index' on 'peptidequants' ('PeptideModSeq' ASC, 'PrecursorCharge' ASC, 'SourceFile' ASC)");
				s.execute("create index if not exists 'PeptideSeq_Peptides_index' on 'peptidequants' ('PeptideSeq' ASC)");

				s.execute(
						"create index if not exists 'PeptideModSeq_PrecursorCharge_SourceFile_Localizations_index' on 'peptidelocalizations' ('PeptideModSeq' ASC, 'PrecursorCharge' ASC, 'SourceFile' ASC)");
				s.execute("create index if not exists 'PeptideSeq_Localizations_index' on 'peptidelocalizations' ('PeptideSeq' ASC)");

				s.execute("create index if not exists 'PeptideModSeq_PrecursorCharge_SourceFile_Scores_index' on 'peptidescores' ('PeptideModSeq' ASC, 'PrecursorCharge' ASC, 'SourceFile' ASC)");
				s.execute("create index if not exists 'PeptideSeq_Scores_index' on 'peptidescores' ('PeptideSeq' ASC)");

				s.execute("create index if not exists 'ProteinGroup_ProteinScores_index' on 'proteinscores' ('ProteinGroup' ASC)");
				s.execute("create index if not exists 'ProteinAccession_ProteinScores_index' on 'proteinscores' ('ProteinAccession' ASC)");

				s.execute("create index if not exists 'PeptideModSeq_PrecursorCharge_SourceFile_Fragments_index' on 'fragmentquants' ('PeptideModSeq' ASC, 'PrecursorCharge' ASC, 'SourceFile' ASC)");
				s.execute("create index if not exists 'PeptideSeq_Fragments_index' on 'fragmentquants' ('PeptideSeq' ASC)");

				s.execute("create index if not exists 'ProteinAccession_PeptideToProtein_index' on 'peptidetoprotein' ('ProteinAccession' ASC)");
				s.execute("create index if not exists 'PeptideSeq_PeptideToProtein_index' on 'peptidetoprotein' ('PeptideSeq' ASC)");

				c.commit();
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	public void close() {
		if (tempFile.exists()&&!tempFile.delete()) {
			Logger.errorLine("Error deleting temp ELIB file!");
		}
	}

}