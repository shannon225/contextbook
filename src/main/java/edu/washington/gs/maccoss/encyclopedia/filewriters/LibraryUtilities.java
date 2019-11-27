package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class LibraryUtilities {
	public static void subsetLibrary(final File saveFile, final float rtMinSec, final float rtMaxSec, final HashSet<String> targets, LibraryFile library)
			throws IOException, SQLException, DataFormatException {
		LibraryFile saveLibrary=new LibraryFile();
		saveLibrary.openFile();
		
		ArrayList<LibraryEntry> toWrite=new ArrayList<>();
		for (LibraryEntry entry : library.getAllEntries(false, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()))) {
			if (rtMinSec<=entry.getRetentionTime()&&rtMaxSec>=entry.getRetentionTime()) {
				if (targets.size()==0||(targets.contains(entry.getPeptideSeq())||targets.contains(entry.getPeptideModSeq()))) {
					toWrite.add(entry);
				}
			}
		}
		Logger.logLine("Found "+toWrite.size()+" peptides from "+targets.size()+" target sequences. Writing to ["+saveFile.getAbsolutePath()+"]...");
		
		saveLibrary.dropIndices();
		saveLibrary.addEntries(toWrite);
		saveLibrary.addProteinsFromEntries(toWrite);
		saveLibrary.createIndices();
		saveLibrary.saveAsFile(saveFile);
		
		saveLibrary.close();
	}
	
	/*
	 * FOR MSPLIT-DIA
	 */
	public static void libraryToMGF(final File saveFile, LibraryFile library)
			throws IOException, SQLException, DataFormatException {
		PrintWriter writer=new PrintWriter(saveFile, "UTF-8");
		
		int count=0;
		for (LibraryEntry entry : library.getAllEntries(false, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()))) {
			count++;
			writer.println("BEGIN IONS");
			writer.println("Title=Scan Number: "+entry.getSpectrumIndex()+" Retention Time: PT"+entry.getScanStartTime()+"S PROTEIN:"+PSMData.accessionsToString(entry.getAccessions()));
			writer.println("CHARGE="+entry.getPrecursorCharge());
			writer.println("PEPMASS="+(entry.getPrecursorMZ()*entry.getPrecursorCharge()-entry.getPrecursorCharge()*MassConstants.protonMass));
			writer.println("SEQ="+entry.getPeptideModSeq().replace("[", "").replace("]", ""));
			writer.println("RTINSECONDS="+entry.getScanStartTime());
			writer.println("SCANS="+Math.round(entry.getScanStartTime())*1000); // assumes scans happen less often than every 1ms

			double[] masses=entry.getMassArray();
			float[] intensities=entry.getIntensityArray();
			for (int i = 0; i < intensities.length; i++) {
				writer.println(masses[i]+" "+intensities[i]);
			}
			
			writer.println("END IONS");
			writer.println();
		}
		writer.close();
		Logger.logLine("Found "+count+" peptides. Writing to ["+saveFile.getAbsolutePath()+"]...");
	}
	
	public static void main(String[] args) throws Exception {
		File libFile=new File("/Volumes/searle_ssd/malaria/novo_yeast/libraries/uniprot_yeast_25jan2019.fasta.z2_nce33.dlib");
		File mgfFile=new File("/Volumes/searle_ssd/malaria/novo_yeast/libraries/uniprot_yeast_25jan2019.fasta.z2_nce33.mgf");
		
		LibraryFile library=new LibraryFile();
		library.openFile(libFile);
		
		libraryToMGF(mgfFile, library);
	}
}
