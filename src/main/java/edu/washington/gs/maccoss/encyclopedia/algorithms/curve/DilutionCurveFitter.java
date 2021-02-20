package edu.washington.gs.maccoss.encyclopedia.algorithms.curve;

import java.awt.Color;
import java.awt.Dimension;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
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
	public static void main2(String[] args) {
		float[] actual = { 4904007.5f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f };
		float[] actual1 = { 6.97E+07f, 4.91E+07f, 3.30E+07f, 1.30E+07f, 6912003f, 3657516.8f, 1903112.4f, 1318391.9f,
				1004388.94f, 850045.44f, 728238.9f, 783248.2f, 896711.94f, 1821167f, 753074.9f, 815791.44f, 551736f,
				461205.78f };
		float[] actual2 = { 7.51E+08f, 6.84E+08f, 4.71E+08f, 1.90E+08f, 7.61E+07f, 3.21E+07f, 1.41E+07f, 6302937.5f,
				2250398f, 822441.75f, 316226.5f, 76445.375f, 334517.88f, 232932.92f, 125391.81f, 350727.62f, 62659.035f,
				582049.1f };
		float[] actual3 = { 1.86E+08f, 1.18E+08f, 7.06E+07f, 2.55E+07f, 1.10E+07f, 5438610f, 2434809f, 1465797.1f,
				124579f, 88526.945f, 70484.29f, 82347.805f, 178821.03f, 73077.97f, 196804.2f, 35501.582f, 90738.1f,
				179157.4f };
		float[] actual4 = { 3187268f, 1503668.9f, 2976294.5f, 1128567f, 265540.84f, 0f, 0f, 153686.38f, 0f, 0f, 0f,
				15378.075f, 0f, 0f, 31187.645f, 0f, 0f, 0f };
		float[] actual5 = { 38458.3f, 9442.072f, 22709.023f, 72533.164f, 72095.195f, 26560.035f, 64255.285f, 38630.027f,
				27547.662f, 11284.466f, 80122.94f, 6767.3096f, 75323.09f, 50457.87f, 23086.375f, 6274.8247f, 16088.956f,
				7947.704f };
		float[] actual6 = { 5255826.5f, 5837896f, 1.08E+08f, 5.33E+07f, 1.84E+07f, 8512274f, 3618325.8f, 1048832f,
				12722.502f, 240561.34f, 90941.555f, 69434.04f, 81829.24f, 69754.37f, 42788.383f, 6405.014f, 93672.79f,
				192462.97f, };
		
		float[] expected = { 1f, 0.68085106f, 0.46666667f, 0.21568628f, 0.1f, 0.04666667f, 0.02156863f, 0.01f,
				0.00466667f, 0.00215686f, 0.001f, 0.00046667f, 0.00021569f, 0.0001f, 4.67E-05f, 2.16E-05f, 0.00001f,
				0.000001f };
		
		expected = new float[] { 1f, 0.1f, 0.01f, 0.001f, 0.000464195f, 0.000215443f, 0.0001f, 4.64E-05f, 1.00E-05f,
				1.00E-06f };
		float[] actual7 = { 1f, 0.065303647f, 0.006984758f, 0.000895221f, 0.000367151f, 0.000143717f, 7.17E-05f,
				3.79E-05f, 0.000141951f, 7.58E-05f };

		TFloatArrayList actualList=new TFloatArrayList(actual7);
		actualList.reverse();
		TFloatArrayList expectedList=new TFloatArrayList(expected);
		expectedList.reverse();
		DilutionFit bestFit=process("PEPTIDE", "PROTEIN", expectedList.toArray(), actualList.toArray()).x;
		ChartPanel panel=graph("PEPTIDE", expectedList.toArray(), actualList.toArray(), bestFit);
		Charter.launchChart(panel, "PEPTIDE");
	}
	
	public static void main(String[] args) throws Exception {
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
		for (FitPeptide fit : fitPeptides) {
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

	private static void addPeptideToAssay(final PrintWriter assayWriter, LibraryEntry entry, float rtInSec, float windowInMin) {
		assayWriter.println(",,(no adduct),"+entry.getPrecursorMZ()+","+entry.getPrecursorCharge()+","+(rtInSec/60f)+","+windowInMin);
	}

	private static float[] incrementDensity(float scanStartTime, float windowInMin, float[] assayDensity) {
		float[] clone=assayDensity.clone();
		int start=Math.round(scanStartTime-windowInMin*60f/2f);
		int stop=Math.round(scanStartTime+windowInMin*60f/2f);
		for (int i = start; i <= stop; i++) {
			if (i<clone.length) {
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
