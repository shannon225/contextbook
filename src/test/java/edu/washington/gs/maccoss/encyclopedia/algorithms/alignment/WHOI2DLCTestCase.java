package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Optional;
import java.util.zip.DataFormatException;

import org.jzy3d.plot3d.builder.Mapper;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter3d;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;

public class WHOI2DLCTestCase {
	private static final PecanSearchParameters PARAMETERS = new PecanSearchParameters(new AminoAcidConstants(),
			FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"),
			false, true, false);
	
	public static void main2(String[] args) throws Exception {

		File reference=new File("/Users/searleb/Downloads/22oct2017_hela_serum_timecourse_narrow_library.elib");
		File twoDLC=new File("/Users/searleb/Downloads/msms.dlib");
		File calibratedFile=new File("/Users/searleb/Downloads/msms_calibrated.dlib");
		
		//File reference=new File("/Volumes/searle_ssd/whoi_bats/190513_1D_BATS_336_DCM_DIA_single_5ug.mzML.elib");
		//File twoDLC=new File("/Volumes/searle_ssd/whoi_bats/190513_1D_BATS_336_DCM_2dDDA.dlib");
		//File calibratedFile=new File("/Volumes/searle_ssd/whoi_bats/calibrated_190513_1D_BATS_336_DCM_2dDDA.dlib");
		
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
		File ref=new File("/Users/searleb/Downloads/22oct2017_hela_serum_timecourse_narrow_library.elib");
		File twoDLC=new File("/Users/searleb/Downloads/msms.dlib");
		//File ref=new File("/Volumes/searle_ssd/whoi_bats/190513_1D_BATS_336_DCM_DIA_single_5ug.mzML.elib");
		//File twoDLC=new File("/Volumes/searle_ssd/whoi_bats/190513_1D_BATS_336_DCM_2dDDA.dlib");
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
		float alpha=Math.min(1.0f, 5000.0f/rts.size());
		XYTraceInterface trace=new XYTrace(rts, GraphType.tinypoint, "Data Used In Fit", new Color(0f, 0f, 1f, alpha), 1.0f);
		Charter.writeAsPDF(new File(twoDLC.getAbsolutePath()+".rt_fit.pdf"), "2D Retention Time (min)", "1D Retention Time (min)", false, trace);
		
		System.out.println("Starting fitting...");
		TwoDimensionalKDE kde=new TwoDimensionalKDE(rts, 1000);

		System.out.println("Starting plotting...");

		final TwoDimensionalKDE finalFilter=kde;
		Mapper mapper=new Mapper() {
			@Override
			public double f(double arg0, double arg1) {
				return finalFilter.f(arg0, arg1);
			}
		};
		Charter3d.plot(mapper, 
				new org.jzy3d.maths.Range(finalFilter.getXRange().getStart(), finalFilter.getXRange().getStop()), 
				new org.jzy3d.maths.Range(finalFilter.getYRange().getStart(), finalFilter.getYRange().getStop()), 
				kde.getResolution()/5);
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
