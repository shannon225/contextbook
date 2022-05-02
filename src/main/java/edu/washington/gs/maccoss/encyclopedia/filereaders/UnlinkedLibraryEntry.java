package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.IOException;
import java.sql.SQLException;
import java.util.HashSet;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class UnlinkedLibraryEntry extends LibraryEntry {
	private static final HashSet<String> UNLOADED_ACCESSIONS=new HashSet<>();
	private final LibraryFile file;
	private HashSet<String> nullableAccessions=null;
	private final boolean isShuffle;
	private final boolean isDecoy;
	private final boolean markAsDecoy;
	private final String originalSequence;
	
	public UnlinkedLibraryEntry(String source, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray, float[] correlationArray, boolean[] quantifiedIonsArray, AminoAcidConstants aaConstants, boolean isShuffle, boolean isDecoy, boolean markAsDecoy, String originalSequence, LibraryFile file) {
		super(source, UNLOADED_ACCESSIONS, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, correlationArray, quantifiedIonsArray, aaConstants);
		this.file=file;
		this.isShuffle=isShuffle;
		this.isDecoy=isDecoy;
		this.markAsDecoy=markAsDecoy;
		this.originalSequence=originalSequence;
	}
	
	public UnlinkedLibraryEntry(String source, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray, AminoAcidConstants aaConstants, boolean isShuffle, boolean isDecoy, boolean markAsDecoy, String originalSequence, LibraryFile file) {
		super(source, UNLOADED_ACCESSIONS, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, aaConstants);
		this.file=file;
		this.isShuffle=isShuffle;
		this.isDecoy=isDecoy;
		this.markAsDecoy=markAsDecoy;
		this.originalSequence=originalSequence;
	}
	
	public void updateAccessions(HashSet<String> accessions) {
		nullableAccessions=new HashSet<String>(accessions);
	}
	
	@Override
	public HashSet<String> getAccessions() {
		if (nullableAccessions==null) {
			try {
				HashSet<String> accessionsFromFile=file.getAccessions(originalSequence);
				if (isShuffle || isDecoy) {
					nullableAccessions = new HashSet<String>();
					for (String accession : accessionsFromFile) {
						if (isShuffle) {
							nullableAccessions.add(SHUFFLE_STRING + accession);
						} else if (isDecoy) {
							nullableAccessions.add(DECOY_STRING + accession);
						}
					}
				} else {
					nullableAccessions = accessionsFromFile;
				}
			}catch (IOException e) {
				throw new EncyclopediaException("Error reading protein accessions from Library File", e);
			} catch (SQLException e) {
				throw new EncyclopediaException("Error reading protein accessions from Library File", e);
			} catch (DataFormatException e) {
				throw new EncyclopediaException("Error reading protein accessions from Library File", e);
			}
		}
		return nullableAccessions;
	}
	
	@Override
	public boolean isDecoy() {
		return markAsDecoy;
	}

	@Override
	protected LibraryEntry updatePeaks(AminoAcidConstants aaConstants, double newPrecursorMz, String newPeptideModSeq, double[] trimmedMasses, float[] trimmedIntensities, float[] trimmedCorrelations, boolean[] trimmedQuantifiedIons, boolean isShuffle, boolean isDecoy, boolean markAsDecoy) {
		return new UnlinkedLibraryEntry(getSource(), getSpectrumIndex(), newPrecursorMz, getPrecursorCharge(), newPeptideModSeq, getCopies(), getRetentionTime(), getScore(), 
				trimmedMasses, trimmedIntensities, trimmedCorrelations, trimmedQuantifiedIons, aaConstants, isShuffle, isDecoy, markAsDecoy, getPeptideSeq(), file);
	}
}
