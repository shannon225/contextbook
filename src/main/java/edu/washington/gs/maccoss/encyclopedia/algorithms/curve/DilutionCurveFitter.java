package edu.washington.gs.maccoss.encyclopedia.algorithms.curve;

import java.awt.Color;
import java.awt.Dimension;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.zip.DataFormatException;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.axis.LogAxis;
import org.jfree.chart.axis.ValueAxis;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearRegression;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class DilutionCurveFitter {
	public static void main(String[] args) {
		
		float[] expected = { 1f, 0.68085106f, 0.46666667f, 0.21568628f, 0.1f, 0.04666667f, 0.02156863f, 0.01f,
				0.00466667f, 0.00215686f, 0.001f, 0.00046667f, 0.00021569f, 0.0001f, 4.67E-05f, 2.16E-05f, 0.00001f};
		
		float[] NLVPMVATVQGQNLK = { 0.000112342f, 6.15181E-05f, 4.9426E-05f, 7.59576E-05f, 0.000100696f, 0.000289449f,
				0.000366852f, 0.001132867f, 0.002648947f, 0.006877976f, 0.014346113f, 0.036599163f, 0.064522792f,
				0.187896813f, 0.36849218f, 0.586212575f, 1f };
		float[] MSSGGGGGDHDHGLSSK = { 0.000111445f, 5.99557E-05f, 5.54356E-05f, 7.87784E-05f, 0.000147337f,
				0.000401697f, 0.000704036f, 0.001509931f, 0.004140635f, 0.008480151f, 0.020829177f, 0.052199693f,
				0.109812594f, 0.194237655f, 0.423290952f, 0.738813441f, 1f };

		float[] VVEQVLR = { 1.22E+08f,	9.25E+07f,	6.25E+07f,	2.59E+07f,	1.17E+07f,	3880509.8f,	1127252.9f,	502679.06f,	199795.1f,	112922.164f,	114000.36f,	143237.23f,	66754.23f,	58291.92f,	33680.6f,	33295.54f,	35985.3f };
		float[] actual=VVEQVLR;//General.reverse(VVEQVLR);

		TFloatArrayList actualList=new TFloatArrayList(actual);
		actualList.reverse();
		TFloatArrayList expectedList=new TFloatArrayList(expected);
		expectedList.reverse();
		DilutionFit bestFit=process("NLVPMVATVQGQNLK", "PROTEIN", expectedList.toArray(), actualList.toArray()).x;
		ChartPanel panel=graph("NLVPMVATVQGQNLK", expectedList.toArray(), actualList.toArray(), bestFit);
		Charter.launchChart(panel, "NLVPMVATVQGQNLK");
	}
	
	public static void main2(String[] args) throws Exception {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		final File outputDirectory=new File("/Users/searleb/Documents/cobbs/2021jan12_cobbs_ln229/curvefitting/");
		final File targetDirectory=new File(outputDirectory, "target");
		final File nontargetDirectory=new File(outputDirectory, "nontarget");
		final File exportLibraryFile=new File(outputDirectory, "target_library.dlib");
		final File libraryAlignmentFile=new File(outputDirectory, "library_rt_alignment.pdf");
		outputDirectory.mkdirs();
		targetDirectory.mkdirs();
		nontargetDirectory.mkdirs();
		
		File dataFile=new File("/Users/searleb/Documents/cobbs/2021jan12_cobbs_ln229/2020dec03_cobbs_cmv_inf_quant.elib.peptides.txt");
		File sampleOrganizationFile=new File("/Users/searleb/Documents/cobbs/2021jan12_cobbs_ln229/sample_organization.csv");
		File libraryFile=new File("/Users/searleb/Documents/cobbs/2021jan12_cobbs_ln229/2020dec03_cobbs_cmv_inf_clib.elib");
		File rtAlignFile=new File("/Users/searleb/Documents/cobbs/2021jan12_cobbs_ln229/2020dec03_cobbs_cmv_curve_dia_0p00_inf.dia.elib");
		int numberOfRTAnchors=10;
		int maxNumberPeptidesPerProtein=3;
		int targetTotalNumberOfPeptides=300; // remember to subtract off anchors (total is 160 peptides)
		float windowInMin=5f; // in minutes!
		final float minCVForAnchors=0.05f;
		final float minCVForBadAnchors=0.75f;
		final int assayMaxDensity=10;
		final String targetAccessionNumberKeyword="HCMV";
		final boolean requireAlignmentRT=true; // turn off for fitting against PRM
		
		final ProgressIndicator progress=new EmptyProgressIndicator();
		
		final HashMap<String, LibraryEntry> libraryEntryByPeptideModSeq=getLibraryData(params, libraryFile);
		final TObjectFloatHashMap<String> knownRTInSecs=new TObjectFloatHashMap<String>();
		ArrayList<XYPoint> rts=new ArrayList<XYPoint>();
		for (Entry<String, LibraryEntry> entry : getLibraryData(params, rtAlignFile).entrySet()) {
			LibraryEntry idealEntry=libraryEntryByPeptideModSeq.get(entry.getKey());
			float alignmentRT = entry.getValue().getScanStartTime();
			knownRTInSecs.put(entry.getKey(), alignmentRT);
			if (idealEntry!=null) {
				XYPoint xy=new XYPoint(idealEntry.getScanStartTime()/60f, alignmentRT/60f);
				rts.add(xy);
			}
		}
		
		if (requireAlignmentRT) {
			for (LibraryEntry entry : new ArrayList<LibraryEntry>(libraryEntryByPeptideModSeq.values())) {
				if (!knownRTInSecs.contains(entry.getPeptideModSeq())) {
					libraryEntryByPeptideModSeq.remove(entry.getPeptideModSeq());
				}
			}
		}
		
		RetentionTimeFilter rtAlignmentFilter=RetentionTimeFilter.getFilter(rts, "Library Retention Time (min)", "Alignment Retention Time (min)");
		rtAlignmentFilter.plot(rts, Optional.of(libraryAlignmentFile));
		final AlignmentWithAnchors rtAlignment=new AlignmentWithAnchors(rtAlignmentFilter, knownRTInSecs);
		
		float minRTInSec=Float.MAX_VALUE;
		float maxRTInSec=-Float.MAX_VALUE;
		for (LibraryEntry entry : libraryEntryByPeptideModSeq.values()) {
			float rtInSec = rtAlignment.getAlignedRTInSec(entry);
			if (rtInSec>maxRTInSec) maxRTInSec=rtInSec;
			if (rtInSec<minRTInSec) minRTInSec=rtInSec;
		}
		Range rtInSecRange=new Range(minRTInSec, maxRTInSec);
		rtInSecRange=new Range(12*60f, 95*60f);
		final ArrayList<Range> subRanges=rtInSecRange.chunkIntoBins(numberOfRTAnchors);
		final String[] bestAnchorPeptideModSeqs=new String[subRanges.size()];
		final float[] bestIntensities=new float[subRanges.size()];
		final float[] bestIntensitiesWithBadCVs=new float[subRanges.size()];
		
		final ArrayList<ScoredObject<String>> expectedConcentrations = getExpectedConcentrationsFromCSV(sampleOrganizationFile);
		final float[] expected = adjustForZeroConcentrations(expectedConcentrations);
		
		final PrintWriter reportWriter=new PrintWriter(new File(outputDirectory, "report.csv"), "UTF-8");
		reportWriter.println("peptide,protein,lod,loq,r2,m,b");

		final ArrayList<FitPeptide> fitPeptides=new ArrayList<FitPeptide>();
		TableParser.parseTSV(dataFile, new TableParserMuscle() {
			public void processRow(Map<String, String> row) {
				String peptide=row.get("Peptide");
				String protein=row.get("Protein");
				
				LibraryEntry entry=libraryEntryByPeptideModSeq.get(peptide);
				
				TFloatArrayList actual=new TFloatArrayList();
				for (ScoredObject<String> scoredObject : expectedConcentrations) {
					String column=scoredObject.y;
					float concentration=Float.parseFloat(row.get(column));
					actual.add(concentration);
				}
				
				float[] actualArray = actual.toArray();
				
				if (protein.indexOf(targetAccessionNumberKeyword)==-1) {
					if (entry!=null) {
						float mean = General.mean(actualArray);
						float cv=General.stdev(actualArray)/mean;
						for (int i = 0; i < bestIntensities.length; i++) {
							float rtInSec = rtAlignment.getAlignedRTInSec(entry);
							if (subRanges.get(i).contains(rtInSec)) {
								if (cv<minCVForAnchors) {
									if (mean>bestIntensities[i]) {
										bestIntensities[i]=mean;
										bestAnchorPeptideModSeqs[i]=peptide;
									}
								} else if (bestIntensities[i]==0.0f&&cv<minCVForBadAnchors) {
									if (mean>bestIntensitiesWithBadCVs[i]&&General.min(actualArray)>0.0f) {
										bestIntensitiesWithBadCVs[i]=mean;
										bestAnchorPeptideModSeqs[i]=peptide;
									}
								}
								break;
							}
						}
					}
					// skip this for curve fitting
					return;
				}
				
				if (entry==null) {
					Logger.errorLine("Found target peptide with no reference RT: "+peptide+", skipping curve fitting!");
					return;
				}

				actualArray=General.divide(actualArray, General.max(actualArray));
				Pair<DilutionFit, Float> pair=process(peptide, protein, expected, actualArray);
				DilutionFit bestFit=pair.x;

				float lod=bestFit.getLOD();
				float loq=bestFit.getLOQ();
				progress.update(peptide+" LOD: "+lod+", LOQ: "+loq);
				
				reportWriter.println(peptide+","+protein+","+lod+","+loq+","+pair.y+","+bestFit.m+","+bestFit.b);
				if (Float.isFinite(loq)&&loq<0) {
					fitPeptides.add(new FitPeptide(peptide, protein, bestFit, expected, actualArray));
				}
			}
			
			public void cleanup() {
			}
		});
		
		reportWriter.flush();
		reportWriter.close();
		
		Collections.sort(fitPeptides);
		
		//////////////
		// BUILD ASSAY
		//////////////

		ArrayList<LibraryEntry> targetEntries=new ArrayList<LibraryEntry>();
		boolean hitMaxDensity=false;
		float[] assayRT=new float[Math.round(rtInSecRange.getStop()+windowInMin*60f)]; // N+W minutes in second increments
		for (int i = 0; i < assayRT.length; i++) {
			assayRT[i]=i/60f;
		}
		float[] assayDensity=new float[assayRT.length];

		final PrintWriter assayWriter=new PrintWriter(new File(outputDirectory, "assay.csv"), "UTF-8");
		assayWriter.println("Compound,Formula,Adduct,m/z,z,RT Time (min),Window (min)");
		
		for (int i = 0; i < bestAnchorPeptideModSeqs.length; i++) {
			if (bestAnchorPeptideModSeqs[i]!=null) {
				LibraryEntry entry=libraryEntryByPeptideModSeq.get(bestAnchorPeptideModSeqs[i]);
				float rtInSec = rtAlignment.getAlignedRTInSec(entry);
				targetEntries.add(entry.updateRetentionTime(rtInSec));

				assayDensity=incrementDensity(rtInSec, windowInMin, assayDensity);
				addPeptideToAssay(assayWriter, entry, rtInSec, windowInMin);
				Logger.logLine("Using "+entry.getPeptideModSeq()+" from "+PSMData.accessionsToString(entry.getAccessions())+" as anchor (rt: "+(rtInSec/60f)+" mins, intensity: "+bestIntensities[i]+" for the RT range from "+(subRanges.get(i).getStart()/60f)+" min to "+(subRanges.get(i).getStop()/60f)+" min");
			} else {
				Logger.logLine("Failed to find good anchor for the RT range from "+(subRanges.get(i).getStart()/60f)+" min to "+(subRanges.get(i).getStop()/60f)+" min");
			}
		}
		
		int count=0;
		HashMap<String, ArrayList<FitPeptide>> targetPeptidesByProtein=new HashMap<String, ArrayList<FitPeptide>>();
		ArrayList<FitPeptide> nontargetedPeptides=new ArrayList<FitPeptide>();
		addpeptides:for (FitPeptide fit : fitPeptides) {
			if (true) break;
			ArrayList<FitPeptide> list=targetPeptidesByProtein.get(fit.proteinKey);
			if (list==null) {
				list=new ArrayList<DilutionCurveFitter.FitPeptide>();
				targetPeptidesByProtein.put(fit.proteinKey, list);
			}

			boolean keep=true;
			if (count<targetTotalNumberOfPeptides) {
				if (list.size()<maxNumberPeptidesPerProtein) {
					LibraryEntry entry=libraryEntryByPeptideModSeq.get(fit.peptideModSeq);

					float rtInSec = rtAlignment.getAlignedRTInSec(entry);
					float[] testDensity=incrementDensity(rtInSec, windowInMin, assayDensity);
					for (int i = 0; i < testDensity.length; i++) {
						if (testDensity[i]>assayMaxDensity) {
							keep=false;
							
							if (!hitMaxDensity) {
								hitMaxDensity=true;
								Logger.logLine("First hit of max density at LOQ: "+fit.bestFit.getLOQ());
							}
							break;
						}
					}
					
					if (keep) {
						assayDensity=testDensity; // update density
						count++;
						Logger.logLine("Adding peptide ("+count+") to assay: "+fit.peptideModSeq+" --> LOQ: "+fit.bestFit.getLOQ()+" from "+fit.proteinKey);
						list.add(fit);
					}
				}
			}
			if (!keep) {
				nontargetedPeptides.add(fit);
			}
		}
		
		count=0;
		int numSingletons=0;
		ArrayList<String> keys=new ArrayList<String>(targetPeptidesByProtein.keySet());
		Collections.sort(keys);
		for (String key : keys) {
			ArrayList<FitPeptide> list=targetPeptidesByProtein.get(key);
			for (FitPeptide fit : list) {
				LibraryEntry entry=libraryEntryByPeptideModSeq.get(fit.peptideModSeq);
				float rtInSec = rtAlignment.getAlignedRTInSec(entry);
				addPeptideToAssay(assayWriter, entry, rtInSec, windowInMin);
				targetEntries.add(entry.updateRetentionTime(rtInSec));

				ChartPanel panel=graph(fit.peptideModSeq, fit.expectedRelativeIntensities, fit.actualRelativeIntensities, fit.bestFit);
				Charter.writeAsPDF(panel.getChart(), new File(targetDirectory, fit.peptideModSeq+".pdf"), new Dimension(300, 300));
				
				count++;
			}
			if (list.size()==1) {
				numSingletons++;
			} if (list.size()==0) {
				targetPeptidesByProtein.remove(key);
			}
		}
		
		for (FitPeptide fit : nontargetedPeptides) {
			ChartPanel panel=graph(fit.peptideModSeq, fit.expectedRelativeIntensities, fit.actualRelativeIntensities, fit.bestFit);
			Charter.writeAsPDF(panel.getChart(), new File(nontargetDirectory, fit.peptideModSeq+".pdf"), new Dimension(300, 300));
		}
		
		assayWriter.flush();
		assayWriter.close();

		LibraryFile exportLibrary=new LibraryFile();
		exportLibrary.openFile();
		exportLibrary.dropIndices();
		exportLibrary.addEntries(targetEntries);
		exportLibrary.addProteinsFromEntries(targetEntries);
		exportLibrary.addMetadata(params.toParameterMap());
		exportLibrary.createIndices();
		exportLibrary.saveAsFile(exportLibraryFile);
		exportLibrary.close();

		XYTrace trace=new XYTrace(assayRT, assayDensity, GraphType.area, "Scheduling density");
		ChartPanel panel=Charter.getChart("Retention Time (min)", "Number of Peptides", true, trace);
		Charter.writeAsPDF(panel.getChart(), new File(outputDirectory, "assay_density.pdf"), new Dimension(600, 300));
		
		Logger.logLine("Finished writing assay for "+targetPeptidesByProtein.size()+" proteins using "+count+" total peptides ("+numSingletons+" single peptide targets)");
		
	}

	protected static void addPeptideToAssay(final PrintWriter assayWriter, LibraryEntry entry, float rtInSec, float windowInMin) {
		assayWriter.println(",,(no adduct),"+entry.getPrecursorMZ()+","+entry.getPrecursorCharge()+","+(rtInSec/60f)+","+windowInMin);
	}

	protected static float[] incrementDensity(float scanStartTime, float windowInMin, float[] assayDensity) {
		float[] clone=assayDensity.clone();
		int start=Math.round(scanStartTime-windowInMin*60f/2f);
		int stop=Math.round(scanStartTime+windowInMin*60f/2f);
		for (int i = start; i <= stop; i++) {
			if (i<clone.length&&i>=0) {
				clone[i]++;
			}
		}
		return clone;
	}
    
	private static HashMap<String, LibraryEntry> getLibraryData(SearchParameters params, File rtAlignFile) throws IOException, SQLException, DataFormatException {
		LibraryFile rtAlignLibrary=new LibraryFile();
		rtAlignLibrary.openFile(rtAlignFile);
		ArrayList<LibraryEntry> entries=rtAlignLibrary.getAllEntries(false, params.getAAConstants());
		HashMap<String, LibraryEntry> rtInSecByPeptideModSeq=new HashMap<String, LibraryEntry>();
		for (LibraryEntry libraryEntry : entries) {
			rtInSecByPeptideModSeq.put(libraryEntry.getPeptideModSeq(), libraryEntry);
		}
		rtAlignLibrary.close();
		return rtInSecByPeptideModSeq;
	}

	private static float[] adjustForZeroConcentrations(final ArrayList<ScoredObject<String>> expectedConcentrations) {
		float minNonZero=Float.MAX_VALUE;
		for (ScoredObject<String> scoredObject : expectedConcentrations) {
			if (scoredObject.getScore()>0&&scoredObject.getScore()<minNonZero) {
				minNonZero=scoredObject.getScore();
			}
		}
		
		Collections.sort(expectedConcentrations);
		TFloatArrayList expectedList=new TFloatArrayList();
		for (ScoredObject<String> scoredObject : expectedConcentrations) {
			float score=scoredObject.getScore();
			if (score==0) score=minNonZero/10f; 
			expectedList.add(score);
		}
		final float[] expected=expectedList.toArray();
		return expected;
	}

	private static ArrayList<ScoredObject<String>> getExpectedConcentrationsFromCSV(File sampleOrganizationFile) {
		final ArrayList<ScoredObject<String>> expectedConcentrations=new ArrayList<ScoredObject<String>>();
		
		System.out.println("Reading "+sampleOrganizationFile.getName()+"...");
		TableParser.parseCSV(sampleOrganizationFile, new TableParserMuscle() {
			public void processRow(Map<String, String> row) {
				String name=row.get("filename");
				float concentration=Float.parseFloat(row.get("concentration"));
				expectedConcentrations.add(new ScoredObject<String>(concentration, name));
			}
			
			public void cleanup() {
			}
		});
		return expectedConcentrations;
	}
	
	public static Pair<DilutionFit, Float> process(String peptide, String protein, float[] expected, float[] actual) {
		TFloatArrayList loggedActual=new TFloatArrayList();
		TFloatArrayList loggedExpected=new TFloatArrayList();

		int startIndex=1; // can't find a crossover below the start index
		float lastZero=-Float.MAX_VALUE;
		float firstNonZero=-Float.MAX_VALUE;
		for (int i = 0; i < actual.length; i++) {
			if (actual[i]>0&&expected[i]>0) {
				float logExpected = Log.log10(expected[i]);
				loggedActual.add(Log.log10(actual[i]));
				loggedExpected.add(logExpected);
				if (firstNonZero==-Float.MAX_VALUE) {
					firstNonZero=logExpected;
				}
			} else {
				if (expected[i]>lastZero) {
					lastZero=expected[i];
					firstNonZero=-Float.MAX_VALUE;
				}
				if (startIndex<loggedActual.size()) {
					startIndex=loggedActual.size()-1; // can't have 0s after the crossover
				}
			}
		}
		if (lastZero>0) {
			// if unspecified, then already set to -maxfloat
			lastZero=Log.log10(lastZero);
		}
		
		if (startIndex==loggedActual.size()) {
			startIndex=0;
		}
		
		float minRSquared=Float.MAX_VALUE;
		DilutionFit fit=null;
		DilutionFit bestFit=null;
		for (int crossOver = startIndex; crossOver < loggedActual.size(); crossOver++) {
			// fit everything below the crossover to a single noise value
			TFloatArrayList noise=new TFloatArrayList();
			for (int j = 0; j <=crossOver; j++) {
				noise.add(loggedActual.get(j));
			}
			float noiseMean = General.mean(noise.toArray());
			
			// fit everything after the crossover to a line
			TFloatArrayList linearX=new TFloatArrayList();
			TFloatArrayList linearY=new TFloatArrayList();
			boolean valuesAboveNoise=false;
			for (int j = crossOver+1; j < loggedActual.size(); j++) {
				linearX.add(loggedExpected.get(j));
				linearY.add(loggedActual.get(j));
				if (loggedActual.get(j)<noiseMean) {
					valuesAboveNoise=true;
				}
			}
			if (valuesAboveNoise) continue;

			// calculate equations
			Pair<Float, Float> equation=LinearRegression.getRegression(linearX.toArray(), linearY.toArray());
			fit=new DilutionFit(noiseMean, General.stdev(noise.toArray()), equation.x, equation.y, lastZero, firstNonZero);
			
			if (false) { // FIXME
				float max=Log.log10(General.max(actual));
				float[] log10Actual = General.subtract(Log.log10(actual), max);
				float[] log10Expected = Log.log10(expected);
				XYTrace values=new XYTrace(log10Expected, log10Actual, GraphType.bigpoint, "Values", Color.black, 4f);
				XYTrace noiseLine=new XYTrace(new float[] {log10Expected[0], log10Expected[expected.length-1]}, new float[] {noiseMean-max, noiseMean-max}, GraphType.dashedline, "Noise", Color.red, 3f);
				XYTrace fitLine=new XYTrace(new float[] {log10Expected[0], log10Expected[expected.length-1]}, new float[] {log10Expected[0]*equation.x+equation.y-max, log10Expected[expected.length-1]*equation.x+equation.y-max}, GraphType.dashedline, "Fit", Color.blue.brighter(), 3f);
				XYTrace pivot=new XYTrace(new float[] {log10Expected[crossOver]}, new float[] {log10Actual[crossOver]}, GraphType.bigpoint, "Pivot", Color.GREEN, 4f);
				ChartPanel panel=Charter.getChart("Expected", "Actual", false, new XYTraceInterface[] {pivot, values, noiseLine, fitLine});
				Charter.launchComponent(panel, "Iteration "+crossOver, new Dimension(300, 300));
			}
			
			if(crossOver>0&&fit.getLOD()<loggedExpected.get(crossOver-1)) {
				// if the point where it hits noiseMean is less than the crossOver point, forcing intercept at noiseMean crossOver point
				equation=LinearRegression.getRegressionWithFixedIntercept(linearX.toArray(), linearY.toArray(), new XYPoint(loggedExpected.get(crossOver), noiseMean));
				fit=new DilutionFit(noiseMean, General.stdev(noise.toArray()), equation.x, equation.y, lastZero, firstNonZero);
			}
			
			// calculate deviation to find the best fit
			float rsquared=0;
			for (int j = 0; j < loggedExpected.size(); j++) {
				float x=loggedExpected.get(j);
				float actualY=loggedActual.get(j);
				float predictedY=fit.getPredicted(x);
				
				float residual=actualY-predictedY;
				rsquared+=residual*residual;
			}
			//System.out.println(crossOver+") "+rsquared+" ("+noise.size()+"/"+linearX.size()+") --> m:"+equation.x+", b:"+equation.y+", lastZero:"+lastZero); // FIXME
			
			// slope has to be at least 0.5
			if (fit.m>=0.5f&&rsquared<minRSquared) {
				minRSquared=rsquared;
				bestFit=fit;
			}
		}
		if (bestFit==null) bestFit=fit; // no good matches, so use last fit
		return new Pair<DilutionFit, Float>(bestFit,minRSquared);
	}
	
	public static ChartPanel graph(String peptide, float[] expected, float[] actual, DilutionFit bestFit) {
		expected=expected.clone();
		actual=actual.clone();
		
		float minNonZeroExpected=Float.MAX_VALUE;
		float minNonZeroActual=Float.MAX_VALUE;
		for (int i = 0; i < actual.length; i++) {
			if (actual[i]>0&&expected[i]>0) {
				if (actual[i]<minNonZeroActual) minNonZeroActual=actual[i];
				if (expected[i]<minNonZeroExpected) minNonZeroExpected=expected[i];
			}
		}
		
		for (int i = 0; i < actual.length; i++) {
			if (actual[i]<=0) actual[i]=minNonZeroActual/10f;
			if (expected[i]<=0) expected[i]=minNonZeroExpected/10f;
		}
		
		float lod=(float)Math.pow(10, bestFit.getLOD());
		float loq=(float)Math.pow(10, bestFit.getLOQ());
		float maxActual=General.max(actual);
		XYTrace lodTrace=new XYTrace(new float[] {lod, lod}, new float[] {minNonZeroActual, maxActual}, GraphType.line, "LOD="+lod, Color.red, 2f);
		XYTrace loqTrace=new XYTrace(new float[] {loq, loq}, new float[] {minNonZeroActual, maxActual}, GraphType.line, "LOQ="+loq, Color.cyan, 2f);
		
		XYTrace actualTrace=new XYTrace(expected, actual, GraphType.bigpoint, peptide, Color.BLACK, 10f);
		XYTrace calculatedTrace=new XYTrace(expected, bestFit.getUnloggedPredicted(expected), GraphType.dashedline, "Calculated", Color.gray, 2f);
		ChartPanel panel=Charter.getChart("Expected", "Actual", true, actualTrace, calculatedTrace, lodTrace, loqTrace);
		
		ValueAxis domain=panel.getChart().getXYPlot().getDomainAxis();
		ValueAxis range=panel.getChart().getXYPlot().getRangeAxis();
		LogAxis newDomain = new LogAxis("Expected");
		newDomain.setLabelFont(domain.getLabelFont());
		newDomain.setTickLabelFont(domain.getTickLabelFont());
		LogAxis newRange =new LogAxis("Actual");
		newRange.setLabelFont(range.getLabelFont());
		newRange.setTickLabelFont(range.getTickLabelFont());

		panel.getChart().getXYPlot().setDomainAxis(newDomain);
		panel.getChart().getXYPlot().setRangeAxis(newRange);
		return panel;
	}

	public static class AlignmentWithAnchors {
		final RetentionTimeFilter rtAlignment;
		final TObjectFloatHashMap<String> knownRTInSecs;
		public AlignmentWithAnchors(RetentionTimeFilter rtAlignment, TObjectFloatHashMap<String> knownRTInSecs) {
			this.rtAlignment = rtAlignment;
			this.knownRTInSecs = knownRTInSecs;
		}

		public float getAlignedRTInSec(LibraryEntry entry) {
			if (knownRTInSecs.contains(entry.getPeptideModSeq())) {
				return knownRTInSecs.get(entry.getPeptideModSeq());
			}
			System.err.println("HAD TO LOOK UP "+entry.getPeptideModSeq());
			return rtAlignment.getYValue(entry.getScanStartTime()/60f)*60f; // deal with sec to min interconversion
		}
	}
	
	public static class FitPeptide implements Comparable<FitPeptide> {
		private final String peptideModSeq;
		private final String proteinKey;
		private final DilutionFit bestFit;
		private float[] expectedRelativeIntensities;
		private float[] actualRelativeIntensities;
		
		public FitPeptide(String peptideModSeq, String proteinKey, DilutionFit bestFit, float[] expectedRelativeIntensities, float[] actualRelativeIntensities) {
			this.peptideModSeq = peptideModSeq;
			this.proteinKey = proteinKey;
			this.bestFit = bestFit;
			this.expectedRelativeIntensities = expectedRelativeIntensities;
			this.actualRelativeIntensities = actualRelativeIntensities;
		}
		
		public int compareTo(FitPeptide o) {
			if (o==null) return 1;
			int c=Float.compare(bestFit.getLOQ(), o.bestFit.getLOQ());
			if (c!=0) return c;

			c=Float.compare(bestFit.getLOD(), o.bestFit.getLOD());
			if (c!=0) return c;
			return peptideModSeq.compareTo(o.peptideModSeq);
		}
	}
	
	public static class DilutionFit {
		private final float noiseValue;
		private final float noiseStdev;
		private final float m;
		private final float b;
		private final float lastZero;
		private final float firstNonZero;
		
		public DilutionFit(float noiseValue, float noiseStdev, float m, float b, float lastZero, float firstNonZero) {
			this.noiseValue = noiseValue;
			this.noiseStdev = noiseStdev;
			this.m = m;
			this.b = b;
			this.lastZero=lastZero;
			this.firstNonZero=firstNonZero;
		}
		
		public float getPredicted(float x) {
			float expectedY=m*x+b;
			if (expectedY<noiseValue) {
				return noiseValue;
			}
			return expectedY;
		}
		
		public float[] getPredicted(float[] xs) {
			float[] expectedYs=new float[xs.length];
			for (int i = 0; i < expectedYs.length; i++) {
				expectedYs[i]=getPredicted(xs[i]);
			}
			return expectedYs;
		}
		
		public float getUnloggedPredicted(float x) {
			float loggedX=Log.log10(x);
			float expectedY=getPredicted(loggedX);
			return (float)Math.pow(10, expectedY);
		}
		
		public float[] getUnloggedPredicted(float[] xs) {
			float[] expectedYs=new float[xs.length];
			for (int i = 0; i < expectedYs.length; i++) {
				expectedYs[i]=getUnloggedPredicted(xs[i]);
			}
			return expectedYs;
		}
		
		public float getLOD() {
			//noiseValue=mx+b
			if (m==0) return Float.POSITIVE_INFINITY;
			return Math.max(lastZero, Math.min(0f, (noiseValue-b)/m));
		}
		public float getLOQ() {
			if (m==0) return Float.POSITIVE_INFINITY;
			float target=noiseValue+3*noiseStdev;
			return Math.max(firstNonZero, Math.min(0f, (target-b)/m));
		}
	}
}
