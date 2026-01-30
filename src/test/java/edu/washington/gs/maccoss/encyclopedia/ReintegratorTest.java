package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ChromatogramLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.ChromatogramExtractor;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearInterpolatedFunction;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LogQuadraticInterpolatedFunction;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;

public class ReintegratorTest {

	public static void main2(String[] args) throws Exception {
		File mainDir=new File("/Users/searle.brian/Documents/temp/mapms/");
		File spectronautCSV=new File(mainDir, "DirectDIA_2024_06_17_125ng_HeLa_1e6_6xpDIA_aurora_01.dia.csv");
		File[] rawFiles=new File[] {new File(mainDir, "2024_06_17_125ng_HeLa_1e6_6xpDIA_aurora_01.dia")};
		//File[] rawFiles=new File[] {new File(mainDir, "2024_06_17_125ng_HeLa_1e6_6xpDIA_aurora_01.dia"), new File(mainDir, "2024_06_17_125ng_HeLa_1e6_6xpDIA_aurora_01_first.dia"), new File(mainDir, "2024_06_17_125ng_HeLa_1e6_6xpDIA_aurora_01_second.dia")};
		File newLibraryFile=new File(mainDir, "2024_06_17_125ng_HeLa_1e6_6xpDIA_aurora_01.elib");
		File prositLibraryFile=new File(mainDir, "BGS_TP_SN19_125ng_1e6_6xpDIA_directDIA_Report.csv_predicted.dlib");
		
		HashMap<String, String> defaultParameters = SearchParameterParser.getDefaultParameters();
		defaultParameters.put(SearchParameters.SMOOTH_INTEGRATIONS, Boolean.FALSE.toString());
		defaultParameters.put(SearchParameters.SUBTRACT_BACKGROUND, Boolean.FALSE.toString());
		final SearchParameters params=SearchParameterParser.parseParameters(defaultParameters);
		
		Reintegrator.reintegrateFromSpectronaut(spectronautCSV, rawFiles, newLibraryFile, Optional.of(prositLibraryFile), params);
	}
	
	public static void main1(String[] args) throws Exception {
		File mainDir=new File("/Users/searle.brian/Documents/temp/mapms/");
		File spectronautCSV=new File(mainDir, "DirectDIA_2024_06_17_125ng_HeLa_1e6_6xpDIA_aurora_01.dia.csv");
		File[] rawFiles=new File[] {new File(mainDir, "2024_06_17_125ng_HeLa_1e6_6xpDIA_aurora_01_first.dia"), new File(mainDir, "2024_06_17_125ng_HeLa_1e6_6xpDIA_aurora_01_second.dia")};
		File prositLibraryFile=new File(mainDir, "BGS_TP_SN19_125ng_1e6_6xpDIA_directDIA_Report.csv_predicted.dlib");
		
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FragmentationType fragType=FragmentationType.HCD;

		LibraryFile prositLibrary=new LibraryFile();
		prositLibrary.openFile(prositLibraryFile);
		ArrayList<LibraryEntry> entries=prositLibrary.getAllEntries(false, parameters.getAAConstants());
		
		StripeFile[] files=new StripeFile[rawFiles.length];
		for (int i=0; i<files.length; i++) {
			files[i]=new StripeFile(true);
			files[i].openFile(rawFiles[i]);
		}
		
		TFloatArrayList deltaTotalIntensityLinear=new TFloatArrayList();
		TFloatArrayList deltaApexIntensityLinear=new TFloatArrayList();
		TFloatArrayList deltaApexRTLinear=new TFloatArrayList();

		TFloatArrayList deltaTotalIntensityLogQuad=new TFloatArrayList();
		TFloatArrayList deltaApexIntensityLogQuad=new TFloatArrayList();
		TFloatArrayList deltaApexRTLogQuad=new TFloatArrayList();
		
		ArrayList<String> peptides=new ArrayList<String>();

		TableParser.parseCSV(spectronautCSV, new TableParserMuscle() {
			
			@Override
			public void processRow(Map<String, String> row) {
				try {
					
					String filename=row.get("R.FileName");
					
					peptides.add(filename);
					if (peptides.size()%1000==0) System.out.print(".");
					if (peptides.size()%10000==0) System.out.print(" ");
					if (peptides.size()%50000==0) System.out.println(" "+peptides.size());
					
//					String quantityString=row.get("FG.Quantity");
//					if (Double.parseDouble(quantityString)<50000) return;
					
					String accessions=row.get("PG.ProteinAccessions");
					HashSet<String> accessionSet=new HashSet<String>(Arrays.asList(accessions.split(";")));
					String intPIMID=row.get("FG.IntMID");
					intPIMID=intPIMID.substring(1, intPIMID.length() - 1);
					String peptideModSeq=PeptideUtils.getCorrectedMasses(intPIMID);
					
					String mzString=row.get("FG.PrecMz");
					double precursorMZ=Double.parseDouble(mzString);
					
					String chargeString=row.get("FG.Charge");
					byte precursorCharge=Byte.parseByte(chargeString);
					
					String rtString=row.get("EG.ApexRT");
					float retentionTime=Float.parseFloat(rtString)*60f; // in sec
					
					String scoreString=row.get("FG.CScore");
					float score=Float.parseFloat(scoreString);
					
					String PEPString=row.get("EG.PEP");
					float PEP=Float.parseFloat(PEPString);
					
					String qvalueString=row.get("EG.Qvalue");
					float qValue=Float.parseFloat(qvalueString);
					
					ArrayList<LibraryEntry> entries=prositLibrary.getEntries(peptideModSeq, precursorCharge, false);
					if (entries.size()==0) return;
					
					AnnotatedLibraryEntry targetEntry=new AnnotatedLibraryEntry(entries.get(0), parameters);
					
					// take the top 6 ions (at least y3)
					FragmentIon[] primaryIonObjects=targetEntry.getMostIntenseAnnotatedIons(3, 3, 0.2f, null);
					
					float[][] totalIntensityLinear=new float[primaryIonObjects.length][];
					float[][] apexIntensityLinear=new float[primaryIonObjects.length][];
					float[][] apexRTLinear=new float[primaryIonObjects.length][];
					float[][] totalIntensityLogQuad=new float[primaryIonObjects.length][];
					float[][] apexIntensityLogQuad=new float[primaryIonObjects.length][];
					float[][] apexRTLogQuad=new float[primaryIonObjects.length][];
					for (int i=0; i<primaryIonObjects.length; i++) {
						totalIntensityLinear[i]=new float[files.length];
						apexIntensityLinear[i]=new float[files.length];
						apexRTLinear[i]=new float[files.length];
						totalIntensityLogQuad[i]=new float[files.length];
						apexIntensityLogQuad[i]=new float[files.length];
						apexRTLogQuad[i]=new float[files.length];
					}
					float halfWidth=15f;	
					
					for (int i=0; i<files.length; i++) {
						StripeFile file=files[i];
						ArrayList<FragmentScan> stripes=file.getStripes(precursorMZ, retentionTime-halfWidth*2, retentionTime+halfWidth*2, false);
	
						ArrayList<Spectrum> downcastedSpectra=FragmentScan.downcastStripeToSpectrum(stripes);
	
						HashMap<FragmentIon, XYTrace> targetFragmentTraceMap=ChromatogramExtractor.extractFragmentChromatograms(parameters.getFragmentTolerance(), primaryIonObjects, downcastedSpectra, null,
								GraphType.boldline, false, false);
						
						for (int j=0; j<primaryIonObjects.length; j++) {
							XYTrace trace=targetFragmentTraceMap.get(primaryIonObjects[j]);
							LinearInterpolatedFunction linearInterpolatedFunction=new LinearInterpolatedFunction(trace.getPoints());
							LogQuadraticInterpolatedFunction logQuadraticInterpolatedFunction=new LogQuadraticInterpolatedFunction(trace.getPoints());
							
							totalIntensityLinear[j][i]=(float)linearInterpolatedFunction.integrate(retentionTime/60f-halfWidth/60f, retentionTime/60f+halfWidth/60f);
							totalIntensityLogQuad[j][i]=(float)logQuadraticInterpolatedFunction.integrate(retentionTime/60f-halfWidth/60f, retentionTime/60f+halfWidth/60f);
							
							XYPoint linearApex=linearInterpolatedFunction.getApex(retentionTime/60f-halfWidth/60f, retentionTime/60f+halfWidth/60f);
							XYPoint logQuadApex=logQuadraticInterpolatedFunction.getApex(retentionTime/60f-halfWidth/60f, retentionTime/60f+halfWidth/60f);
							
							apexIntensityLinear[j][i]=(float)linearApex.y;
							apexIntensityLogQuad[j][i]=(float)logQuadApex.y;
							
							apexRTLinear[j][i]=(float)linearApex.x;
							apexRTLogQuad[j][i]=(float)logQuadApex.x;
						}
					}
					
					for (int j=0; j<primaryIonObjects.length; j++) {
						deltaTotalIntensityLinear.add(Math.abs(totalIntensityLinear[j][0]-totalIntensityLinear[j][1])/Math.max(totalIntensityLinear[j][0], totalIntensityLinear[j][1]));
						deltaApexIntensityLinear.add(Math.abs(apexIntensityLinear[j][0]-apexIntensityLinear[j][1])/Math.max(apexIntensityLinear[j][0], apexIntensityLinear[j][1]));
						deltaApexRTLinear.add(Math.abs(retentionTime/60f-apexRTLinear[j][0]));
						deltaApexRTLinear.add(Math.abs(retentionTime/60f-apexRTLinear[j][1]));

						deltaTotalIntensityLogQuad.add(Math.abs(totalIntensityLogQuad[j][0]-totalIntensityLogQuad[j][1])/Math.max(totalIntensityLogQuad[j][0], totalIntensityLogQuad[j][1]));
						deltaApexIntensityLogQuad.add(Math.abs(apexIntensityLogQuad[j][0]-apexIntensityLogQuad[j][1])/Math.max(apexIntensityLogQuad[j][0], apexIntensityLogQuad[j][1]));
						deltaApexRTLogQuad.add(Math.abs(retentionTime/60f-apexRTLogQuad[j][0]));
						deltaApexRTLogQuad.add(Math.abs(retentionTime/60f-apexRTLogQuad[j][1]));
					}
					
				} catch (Exception e) {
					Logger.errorException(e);
				}
			}
			
			@Override
			public void cleanup() {
			}
		});
		System.out.println();
		System.out.println("deltaTotalIntensityLinear: "+QuickMedian.median(deltaTotalIntensityLinear.toArray()));
		System.out.println("deltaApexIntensityLinear: "+QuickMedian.median(deltaApexIntensityLinear.toArray()));
		System.out.println("deltaApexRTLinear: "+QuickMedian.median(deltaApexRTLinear.toArray()));
		System.out.println();
		System.out.println("deltaTotalIntensityLogQuad: "+QuickMedian.median(deltaTotalIntensityLogQuad.toArray()));
		System.out.println("deltaApexIntensityLogQuad: "+QuickMedian.median(deltaApexIntensityLogQuad.toArray()));
		System.out.println("deltaApexRTLogQuad: "+QuickMedian.median(deltaApexRTLogQuad.toArray()));
	}

	public static void main(String[] args) throws Exception {
		File mainDir=new File("/Users/searle.brian/Documents/temp/mapms_GPFDIA/raws/");
		File searchResults=new File(mainDir, "2024_06_17_125ng_HeLa_6xpGPFDIA_aurora_0combined.dia.encyclopedia2.txt.elib");
		File[] rawFiles=new File[] {new File(mainDir, "2024_06_17_125ng_HeLa_6xpGPFDIA_aurora_0combined_first.dia"), new File(mainDir, "2024_06_17_125ng_HeLa_6xpGPFDIA_aurora_0combined_second.dia")};
		
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		FragmentationType fragType=FragmentationType.HCD;

		LibraryFile searchResultFile=new LibraryFile();
		searchResultFile.openFile(searchResults);
		ArrayList<LibraryEntry> entries=searchResultFile.getAllEntries(false, parameters.getAAConstants());
		
		StripeFile[] files=new StripeFile[rawFiles.length];
		for (int i=0; i<files.length; i++) {
			files[i]=new StripeFile(true);
			files[i].openFile(rawFiles[i]);
		}
		
		TIntObjectHashMap<TFloatArrayList> deltaTotalIntensityLinear=new TIntObjectHashMap<TFloatArrayList>();
		TIntObjectHashMap<TFloatArrayList> deltaApexIntensityLinear=new TIntObjectHashMap<TFloatArrayList>();
		TIntObjectHashMap<TFloatArrayList> deltaApexRTLinear=new TIntObjectHashMap<TFloatArrayList>();

		TIntObjectHashMap<TFloatArrayList> deltaTotalIntensityLogQuad=new TIntObjectHashMap<TFloatArrayList>();
		TIntObjectHashMap<TFloatArrayList> deltaApexIntensityLogQuad=new TIntObjectHashMap<TFloatArrayList>();
		TIntObjectHashMap<TFloatArrayList> deltaApexRTLogQuad=new TIntObjectHashMap<TFloatArrayList>();
		
		Statement s=searchResultFile.getConnection().createStatement();
		for (LibraryEntry entry : entries) {
			AnnotatedLibraryEntry targetEntry=new AnnotatedLibraryEntry(entry, parameters);
			double precursorMZ=targetEntry.getPrecursorMZ();
			
			String sql="SELECT RTInSecondsStart, RTInSecondsCenter, RTInSecondsStop FROM peptidequants WHERE peptidemodseq is \""+entry.getPeptideModSeq()+"\"";
			ResultSet rs=s.executeQuery(sql);
			float rtInSecondsStart=rs.getFloat(1);
			float retentionTime=rs.getFloat(2);
			float rtInSecondsStop=rs.getFloat(3);
			float peakWidthInSec=rtInSecondsStop-rtInSecondsStart;
			if (peakWidthInSec<20) continue;
			
			// take the top 6 ions (at least y3)
			FragmentIon[] primaryIonObjects=targetEntry.getMostIntenseAnnotatedIons(6, 3, 0.2f, entry.getQuantifiedIonsArray());
			//FragmentIon[] primaryIonObjects=targetEntry.getMostCorrelatedAnnotatedIons(6, 3, TransitionRefiner.quantitativeCorrelationThreshold);

			int[][] numberOfPoints=new int[primaryIonObjects.length][];
			float[][] totalIntensityLinear=new float[primaryIonObjects.length][];
			float[][] apexIntensityLinear=new float[primaryIonObjects.length][];
			float[][] apexRTLinear=new float[primaryIonObjects.length][];
			float[][] totalIntensityLogQuad=new float[primaryIonObjects.length][];
			float[][] apexIntensityLogQuad=new float[primaryIonObjects.length][];
			float[][] apexRTLogQuad=new float[primaryIonObjects.length][];
			for (int i=0; i<primaryIonObjects.length; i++) {
				numberOfPoints[i]=new int[files.length];
				totalIntensityLinear[i]=new float[files.length];
				apexIntensityLinear[i]=new float[files.length];
				apexRTLinear[i]=new float[files.length];
				totalIntensityLogQuad[i]=new float[files.length];
				apexIntensityLogQuad[i]=new float[files.length];
				apexRTLogQuad[i]=new float[files.length];
			}	
			//System.out.println(rtInSecondsStart/60f+" --> "+rtInSecondsStop/60f);
			for (int i=0; i<files.length; i++) {
				StripeFile file=files[i];
				ArrayList<FragmentScan> stripes=file.getStripes(precursorMZ, rtInSecondsStart-15, rtInSecondsStop+15, false);

				ArrayList<Spectrum> downcastedSpectra=FragmentScan.downcastStripeToSpectrum(stripes);

				HashMap<FragmentIon, XYTrace> targetFragmentTraceMap=ChromatogramExtractor.extractFragmentChromatograms(parameters.getFragmentTolerance(), primaryIonObjects, downcastedSpectra, null,
						GraphType.boldline, false, false);
				
				for (int j=0; j<primaryIonObjects.length; j++) {
					XYTrace trace=targetFragmentTraceMap.get(primaryIonObjects[j]);
					//trace=SkylineSGFilter.adjustableSavitzkyGolaySmooth(trace);
					
					Range rtRange=new Range(rtInSecondsStart/60f, rtInSecondsStop/60f);
					XYPoint apex=trace.getMaxXYInRange(rtRange);

					int numPoints=0; // points in range above 5% apex
					for (XYPoint point : trace.getPoints()) {
						if (rtRange.contains(point.x)) {
							if (point.y/apex.y>0.05) {
								numPoints++;
							}
						}
					} 
					numberOfPoints[j][i]=numPoints;
					
					LinearInterpolatedFunction linearInterpolatedFunction=new LinearInterpolatedFunction(trace.getPoints());
					LogQuadraticInterpolatedFunction logQuadraticInterpolatedFunction=new LogQuadraticInterpolatedFunction(trace.getPoints(), 1000.0);
					
					totalIntensityLinear[j][i]=(float)linearInterpolatedFunction.integrate(rtRange.getStart(), rtRange.getStop());
					totalIntensityLogQuad[j][i]=(float)logQuadraticInterpolatedFunction.integrate(rtRange.getStart(), rtRange.getStop());
					
					XYPoint linearApex=linearInterpolatedFunction.getApex(rtRange.getStart(), rtRange.getStop());
					XYPoint logQuadApex=logQuadraticInterpolatedFunction.getApex(rtRange.getStart(), rtRange.getStop());
					
					apexIntensityLinear[j][i]=(float)linearApex.y;
					apexIntensityLogQuad[j][i]=(float)logQuadApex.y;
					
					apexRTLinear[j][i]=(float)linearApex.x;
					apexRTLogQuad[j][i]=(float)logQuadApex.x;
					
//					if (numPoints==12&&Math.abs(retentionTime/60f-linearInterpolatedFunction.getApex(rtRange.getStart(), rtRange.getStop()).x)>0.15) {
//						System.out.println(entry.getPeptideModSeq());
//						System.out.println(retentionTime/60f+"\t"+-1);
//						for (XYPoint point : trace.getPoints()) {
//							System.out.println(point.x+"\t"+point.y);
//						}
//						System.exit(1);
//					}
				}
			}
			
			for (int j=0; j<primaryIonObjects.length; j++) {
				int numPoints=Math.min(numberOfPoints[j][0], numberOfPoints[j][1]);
				
				TFloatArrayList deltaTotalIntensityLinearList=deltaTotalIntensityLinear.get(numPoints);
				if (deltaTotalIntensityLinearList==null) {
					deltaTotalIntensityLinearList=new TFloatArrayList();
					deltaTotalIntensityLinear.put(numPoints, deltaTotalIntensityLinearList);
				}
				
				TFloatArrayList deltaApexIntensityLinearList=deltaApexIntensityLinear.get(numPoints);
				if (deltaApexIntensityLinearList==null) {
					deltaApexIntensityLinearList=new TFloatArrayList();
					deltaApexIntensityLinear.put(numPoints, deltaApexIntensityLinearList);
				}
				
				TFloatArrayList deltaApexRTLinearList=deltaApexRTLinear.get(numPoints);
				if (deltaApexRTLinearList==null) {
					deltaApexRTLinearList=new TFloatArrayList();
					deltaApexRTLinear.put(numPoints, deltaApexRTLinearList);
				}

				
				TFloatArrayList deltaTotalIntensityLogQuadList=deltaTotalIntensityLogQuad.get(numPoints);
				if (deltaTotalIntensityLogQuadList==null) {
					deltaTotalIntensityLogQuadList=new TFloatArrayList();
					deltaTotalIntensityLogQuad.put(numPoints, deltaTotalIntensityLogQuadList);
				}
				
				TFloatArrayList deltaApexIntensityLogQuadList=deltaApexIntensityLogQuad.get(numPoints);
				if (deltaApexIntensityLogQuadList==null) {
					deltaApexIntensityLogQuadList=new TFloatArrayList();
					deltaApexIntensityLogQuad.put(numPoints, deltaApexIntensityLogQuadList);
				}
				
				TFloatArrayList deltaApexRTLogQuadList=deltaApexRTLogQuad.get(numPoints);
				if (deltaApexRTLogQuadList==null) {
					deltaApexRTLogQuadList=new TFloatArrayList();
					deltaApexRTLogQuad.put(numPoints, deltaApexRTLogQuadList);
				}
				
				deltaTotalIntensityLinearList.add(Math.abs(totalIntensityLinear[j][0]-totalIntensityLinear[j][1])/Math.max(totalIntensityLinear[j][0], totalIntensityLinear[j][1]));
				deltaApexIntensityLinearList.add(Math.abs(apexIntensityLinear[j][0]-apexIntensityLinear[j][1])/Math.max(apexIntensityLinear[j][0], apexIntensityLinear[j][1]));
				deltaApexRTLinearList.add(100.0f*Math.abs(retentionTime/60f-apexRTLinear[j][0])/peakWidthInSec);
				deltaApexRTLinearList.add(100.0f*Math.abs(retentionTime/60f-apexRTLinear[j][1])/peakWidthInSec);

				deltaTotalIntensityLogQuadList.add(Math.abs(totalIntensityLogQuad[j][0]-totalIntensityLogQuad[j][1])/Math.max(totalIntensityLogQuad[j][0], totalIntensityLogQuad[j][1]));
				deltaApexIntensityLogQuadList.add(Math.abs(apexIntensityLogQuad[j][0]-apexIntensityLogQuad[j][1])/Math.max(apexIntensityLogQuad[j][0], apexIntensityLogQuad[j][1]));
				deltaApexRTLogQuadList.add(100.0f*Math.abs(retentionTime/60f-apexRTLogQuad[j][0])/peakWidthInSec);
				deltaApexRTLogQuadList.add(100.0f*Math.abs(retentionTime/60f-apexRTLogQuad[j][1])/peakWidthInSec);
			}
		}
		
		System.out.println();
		System.out.println("N: ");
		int minN=3;
		int maxN=12;
		for (int i=minN; i<=maxN; i++) {
			System.out.println(i+"\t"+deltaTotalIntensityLinear.get(i-1).size());
		}
		
		System.out.println();
		System.out.println("deltaTotalIntensityLinear\t5p\t25p\t50p\t75p\t95p");
		for (int i=minN; i<=maxN; i++) {
			TIntObjectHashMap<TFloatArrayList> source=deltaTotalIntensityLinear;
			System.out.println(i+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.05f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.25f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.5f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.75f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.95f));
		}
		
		System.out.println("deltaApexIntensityLinear\t5p\t25p\t50p\t75p\t95p");
		for (int i=minN; i<=maxN; i++) {
			TIntObjectHashMap<TFloatArrayList> source=deltaApexIntensityLinear;
			System.out.println(i+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.05f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.25f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.5f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.75f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.95f));
		}
		System.out.println("percentDeltaApexRTLinear\t5p\t25p\t50p\t75p\t95p");
		for (int i=minN; i<=maxN; i++) {
			TIntObjectHashMap<TFloatArrayList> source=deltaApexRTLinear;
			System.out.println(i+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.05f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.25f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.5f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.75f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.95f));
		}
		System.out.println();
		System.out.println("deltaTotalIntensityLogQuad\t5p\t25p\t50p\t75p\t95p");
		for (int i=minN; i<=maxN; i++) {
			TIntObjectHashMap<TFloatArrayList> source=deltaTotalIntensityLogQuad;
			System.out.println(i+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.05f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.25f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.5f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.75f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.95f));
		}
		System.out.println("deltaApexIntensityLogQuad\t5p\t25p\t50p\t75p\t95p");
		for (int i=minN; i<=maxN; i++) {
			TIntObjectHashMap<TFloatArrayList> source=deltaApexIntensityLogQuad;
			System.out.println(i+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.05f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.25f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.5f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.75f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.95f));
		}
		System.out.println("percentDeltaApexRTLogQuad\t5p\t25p\t50p\t75p\t95p");
		for (int i=minN; i<=maxN; i++) {
			TIntObjectHashMap<TFloatArrayList> source=deltaApexRTLogQuad;
			System.out.println(i+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.05f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.25f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.5f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.75f)+"\t"+QuickMedian.select(source.get(i-1).toArray(), 0.95f));
		}
	}
}
