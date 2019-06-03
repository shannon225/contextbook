package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.zip.DataFormatException;

import com.fasterxml.jackson.databind.JsonMappingException.Reference;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface.AlignmentDataPoint;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.QuantitativeSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.LibraryUtilities;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter3d;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;

public class WHOI2DLCTestCase {
	private static final PecanSearchParameters PARAMETERS = new PecanSearchParameters(new AminoAcidConstants(),
			FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"),
			false, true, false);
	
	public static void main2(String[] args) throws Exception {

		//File reference=new File("/Users/searleb/Downloads/22oct2017_hela_serum_timecourse_narrow_library.elib");
		//File twoDLC=new File("/Users/searleb/Downloads/msms.dlib");
		//File calibratedFile=new File("/Users/searleb/Downloads/msms_calibrated.dlib");
		
		File reference=new File("/Volumes/searle_ssd/whoi_bats/190513_1D_BATS_336_DCM_DIA_single_5ug.mzML.elib");
		File twoDLC=new File("/Volumes/searle_ssd/whoi_bats/190513_1D_BATS_336_DCM_2dDDA.dlib");
		File calibratedFile=new File("/Volumes/searle_ssd/whoi_bats/calibrated_190513_1D_BATS_336_DCM_2dDDA.dlib");
		
		LibraryFile referenceLibrary=new LibraryFile();
		referenceLibrary.openFile(reference);
		ArrayList<LibraryEntry> entries=referenceLibrary.getAllEntries(false, PARAMETERS.getAAConstants());
		TObjectFloatHashMap<String> referenceRTInSecs=new TObjectFloatHashMap<>();
		for (LibraryEntry entry : entries) {
			referenceRTInSecs.put(entry.getPeptideModSeq(), entry.getRetentionTime()/60f);
		}

		LibraryFile library=new LibraryFile();
		library.openFile(twoDLC);
		ArrayList<LibraryEntry> allEntries=library.getAllEntries(false, PARAMETERS.getAAConstants());
		library.close();

		ArrayList<LibraryEntry> calibratedEntries=new ArrayList<>();
		float[] timeBoundaries=new float[] {0f, 2f, 32f, 62f, 92f, 122f, 152f, 182f, 212f, 242f, 272f, 302f, 332f, 362f, 392f, 1000f};
		for (int i=1; i<timeBoundaries.length; i++) {
			String key=Math.round(timeBoundaries[i-1])+" to "+Math.round(timeBoundaries[i]);
			System.out.println("writing "+key);
			
			float rtMinSec=timeBoundaries[i-1]*60f;
			float rtMaxSec=timeBoundaries[i]*60f;
			ArrayList<LibraryEntry> inThisWindow=new ArrayList<>();
			for (LibraryEntry entry : allEntries) {
				if (rtMinSec<=entry.getRetentionTime()&&rtMaxSec>=entry.getRetentionTime()) {
					inThisWindow.add(entry);
				}
			}

			TObjectFloatHashMap<String> subsetRTInSecs=new TObjectFloatHashMap<>();
			for (LibraryEntry entry : inThisWindow) {
				subsetRTInSecs.put(entry.getPeptideModSeq(), entry.getRetentionTime()/60f);
			}

			ArrayList<XYPoint> points = ReferencePeakIntegrator.getMatchingPoints(referenceRTInSecs, subsetRTInSecs);

			RetentionTimeAlignmentInterface alignment=RetentionTimeFilter.getFilter(points, "reference", key, 2000);
			File saveFileSeed=new File(twoDLC.getParentFile(), "msms_"+Math.round(timeBoundaries[i-1])+"to"+Math.round(timeBoundaries[i])+".dlib");
			alignment.plot(points, Optional.ofNullable(saveFileSeed)); // save PDF to file

			for (LibraryEntry entry : inThisWindow) {
				LibraryEntry rtCorrectedEntry=entry.updateRetentionTime(60f*alignment.getXValue(entry.getRetentionTime()/60f));
				calibratedEntries.add(rtCorrectedEntry);
			}
		}

		LibraryFile calibrated=new LibraryFile();
		calibrated.openFile();
		Logger.logLine("Found "+calibratedEntries.size()+" peptides. Writing to ["+calibratedFile.getAbsolutePath()+"]...");
		
		calibrated.dropIndices();
		calibrated.addEntries(calibratedEntries);
		calibrated.addProteinsFromEntries(calibratedEntries);
		calibrated.createIndices();
		calibrated.saveAsFile(calibratedFile);
		calibrated.close();
	}

	public static void main(String[] args) throws Exception {
		//File ref=new File("/Users/searleb/Downloads/22oct2017_hela_serum_timecourse_narrow_library.elib");
		//File twoDLC=new File("/Users/searleb/Downloads/msms.dlib");
		File ref=new File("/Volumes/searle_ssd/whoi_bats/190513_1D_BATS_336_DCM_DIA_single_5ug.mzML.elib");
		File twoDLC=new File("/Volumes/searle_ssd/whoi_bats/190513_1D_BATS_336_DCM_2dDDA.dlib");
		//File twoDLC=new File("/Volumes/searle_ssd/whoi_bats/calibrated_190513_1D_BATS_336_DCM_2dDDA.dlib");
		
		LibraryFile refLib=new LibraryFile();
		refLib.openFile(ref);
		LibraryFile twoDLCLib=new LibraryFile();
		twoDLCLib.openFile(twoDLC);
		
		TObjectFloatHashMap<String> refRTsBySeq=getRTs(refLib);
		TObjectFloatHashMap<String> tarRTsBySeq=getRTs(twoDLCLib);
		
		ArrayList<XYPoint> rts=new ArrayList<>();
		
		tarRTsBySeq.forEachEntry(new TObjectFloatProcedure<String>() {
			@Override
			public boolean execute(String arg0, float arg1) {
				float rt=refRTsBySeq.get(arg0);
				if (rt!=0.0f) {
					rts.add(new XYPoint(arg1/60f, rt/60f));
				}
				return true;
			}
		});
		
//		for (XYPoint xyPoint : rts) {
//			System.out.println(xyPoint.toString());
//		}
//		if (true) System.exit(1);

		Collections.sort(rts);
		System.out.println("Starting fitting...");
		TwoDimensionalKDE kde=new TwoDimensionalKDE(rts, 1000);

		System.out.println("Starting plotting...");
		Charter3d.plot(kde, kde.getXRange(), kde.getYRange(), kde.getResolution()/5);
	}

	private static TObjectFloatHashMap<String> getRTs(LibraryFile refLib) throws IOException, SQLException, DataFormatException {
		ArrayList<LibraryEntry> entries=refLib.getAllEntries(false, PARAMETERS.getAAConstants());
		TObjectFloatHashMap<String> rtsBySeq=new TObjectFloatHashMap<>();
		for (LibraryEntry entry : entries) {
			rtsBySeq.put(entry.getPeptideModSeq(), entry.getRetentionTime());
		}
		return rtsBySeq;
	}

}
