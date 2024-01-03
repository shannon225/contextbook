package edu.washington.gs.maccoss.encyclopedia.algorithms.curve;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.sql.SQLException;
import java.text.AttributedString;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.TreeMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.zip.DataFormatException;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.axis.LogAxis;
import org.jfree.chart.axis.LogTick;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.DeviationRenderer;
import org.jfree.chart.renderer.xy.XYDifferenceRenderer;
import org.jfree.data.general.Dataset;
import org.jfree.data.statistics.BoxAndWhiskerItem;
import org.jfree.data.xy.DefaultIntervalXYDataset;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.TextAnchor;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Boxplotter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.ExtendedLogAxis;
import edu.washington.gs.maccoss.encyclopedia.gui.general.NumberBoxAndWhiskerXYDataset;
import edu.washington.gs.maccoss.encyclopedia.gui.general.XYGraphingTrace;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Boxplotter.XYBoxPlotterRenderer;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
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
import gnu.trove.procedure.TObjectFloatProcedure;

public class DilutionCurveFitter {
	public static void main3(String[] args) {
		
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
		DilutionFit bestFit=process("NLVPMVATVQGQNLK", "PROTEIN", expectedList.toArray(), actualList.toArray(), General.max(actual), true).x;
		ChartPanel panel=graph("NLVPMVATVQGQNLK", expectedList.toArray(), actualList.toArray(), bestFit, Optional.empty());
		Charter.launchChart(panel, "NLVPMVATVQGQNLK");
	}
	
	public static void main(String[] args) throws Exception {
		final File outputDirectory=new File("/Users/searleb/Documents/cobbs/2021jan12_cobbs_ln229/figures/prms/curvefitting/");
		final File targetDirectory=new File(outputDirectory, "target");
		outputDirectory.mkdirs();
		targetDirectory.mkdirs();
		
		File dataFile=new File("/Users/searleb/Documents/cobbs/2021jan12_cobbs_ln229/figures/prms/final_119_and_248_PRM_peptide_quant_report.csv");
		File sampleOrganizationFile=new File("/Users/searleb/Documents/cobbs/2021jan12_cobbs_ln229/figures/prms/prm_sample_organization.csv");

		
		Pair<ArrayList<ScoredObject<String>>, Map<String, TObjectFloatHashMap<String>>> concentrationPair=getExpectedConcentrationsFromCSV(sampleOrganizationFile);
		final ArrayList<ScoredObject<String>> expectedConcentrations=concentrationPair.x;
		final Map<String, TObjectFloatHashMap<String>> unknowns=concentrationPair.y;
		final float[] expected = adjustForZeroConcentrations(expectedConcentrations);

		final ArrayList<FitPeptide> fitPeptides=fitCurves(outputDirectory, dataFile, expectedConcentrations, expected, new DilutionCurveFitting8HzWideParameters(), true);
		final HashMap<String, Map<String, TObjectFloatHashMap<String>>> unknownData=extractUnknowns(dataFile, unknowns, "HCMV");
		
		for (FitPeptide fit : fitPeptides) {
			Map<String, TObjectFloatHashMap<String>> data=unknownData.get(fit.peptideModSeq);
			ChartPanel panel=graph(fit.peptideModSeq, fit.expectedRelativeIntensities, fit.actualRelativeIntensities, fit.bestFit, Optional.ofNullable(data));
			Charter.writeAsPDF(panel.getChart(), new File(targetDirectory, fit.peptideModSeq+".pdf"), new Dimension(400, 300));
		}
	}
	
	public static void main2(String[] args) throws Exception {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		final File outputDirectory=new File("/Users/searleb/Documents/cobbs/2021jan12_cobbs_ln229/test/curvefitting_wide_testing/");
		File dataFile=new File("/Users/searleb/Documents/cobbs/2021jan12_cobbs_ln229/test/2020dec03_cobbs_cmv_inf_quant.elib.peptides.txt");
		File sampleOrganizationFile=new File("/Users/searleb/Documents/cobbs/2021jan12_cobbs_ln229/test/sample_organization.csv");
		File libraryFile=new File("/Users/searleb/Documents/cobbs/2021jan12_cobbs_ln229/test/2020dec03_cobbs_cmv_inf_clib.elib");
		File rtAlignFile=new File("/Users/searleb/Documents/cobbs/2021jan12_cobbs_ln229/test/2020dec03_cobbs_cmv_curve_dia_0p00_inf.dia.elib");
		
		//AbstractDilutionCurveFittingParameters fittingParams=new DilutionCurveFitting4HzDeepParameters();
		AbstractDilutionCurveFittingParameters fittingParams=new DilutionCurveFitting8HzWideParameters();
		
		generateAssayFromCurves(params, outputDirectory, dataFile, sampleOrganizationFile, libraryFile, rtAlignFile, fittingParams);
	}

	public static void generateAssayFromCurves(SearchParameters params, final File outputDirectory, File dataFile,
			File sampleOrganizationFile, File libraryFile, File rtAlignFile, AbstractDilutionCurveFittingParameters fittingParams)
			throws IOException, SQLException, DataFormatException, FileNotFoundException, UnsupportedEncodingException {
		final File targetDirectory=new File(outputDirectory, "target");
		final File nontargetDirectory=new File(outputDirectory, "nontarget");
		final File exportLibraryFile=new File(outputDirectory, "target_library.dlib");
		final File libraryAlignmentFile=new File(outputDirectory, "library_rt_alignment.pdf");
		outputDirectory.mkdirs();
		targetDirectory.mkdirs();
		nontargetDirectory.mkdirs();
		
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
		if (fittingParams.isRequireAlignmentRT()) {
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
		//rtInSecRange=new Range(12*60f, 95*60f);
		final ArrayList<Range> subRanges=rtInSecRange.chunkIntoBins(fittingParams.getNumberOfRTAnchors());
		
		Pair<ArrayList<ScoredObject<String>>, Map<String, TObjectFloatHashMap<String>>> concentrationPair= getExpectedConcentrationsFromCSV(sampleOrganizationFile);
		final ArrayList<ScoredObject<String>> expectedConcentrations=concentrationPair.x;
		
		final float[] expected = adjustForZeroConcentrations(expectedConcentrations);
		
		Pair<String[], float[]> anchorData=extractAnchorPeptides(dataFile, fittingParams,
				libraryEntryByPeptideModSeq, rtAlignment, subRanges, expectedConcentrations);

		final ArrayList<FitPeptide> fitPeptides=fitCurves(outputDirectory, dataFile, 
				expectedConcentrations, expected, fittingParams, fittingParams.isUseLineNoise());

		final String[] bestAnchorPeptideModSeqs=anchorData.x;
		final float[] bestIntensities=anchorData.y;

		ArrayList<LibraryEntry> targetEntries = scheduleAssay(outputDirectory, targetDirectory, nontargetDirectory,
				fittingParams, libraryEntryByPeptideModSeq, rtAlignment, rtInSecRange, subRanges,
				bestAnchorPeptideModSeqs, bestIntensities, fitPeptides);

		writeLibraryEntries(params, exportLibraryFile, targetEntries);
	}
	
	private static Pair<String[], float[]> extractAnchorPeptides(File dataFile, AbstractDilutionCurveFittingParameters fittingParams,
			final HashMap<String, LibraryEntry> libraryEntryByPeptideModSeq, final AlignmentWithAnchors rtAlignment,
			final ArrayList<Range> subRanges, final ArrayList<ScoredObject<String>> expectedConcentrations) throws FileNotFoundException, UnsupportedEncodingException {
		
		final String[] bestAnchorPeptideModSeqs=new String[subRanges.size()];
		final float[] bestIntensities=new float[subRanges.size()];
		final float[] bestIntensitiesWithBadCVs=new float[subRanges.size()];
		Pair<String[], float[]> data=new Pair<String[], float[]>(bestAnchorPeptideModSeqs, bestIntensities);

		TableParserMuscle muscle = new TableParserMuscle() {
			public void processRow(Map<String, String> row) {
				String peptide=row.get("Peptide"); // from EncyclopeDIA
				if (peptide==null) peptide=row.get("Peptide Modified Sequence"); // from Skyline
				String protein=row.get("Protein"); // from EncyclopeDIA
				if (protein==null) protein=row.get("Protein Name"); // from Skyline
				
				LibraryEntry entry=libraryEntryByPeptideModSeq.get(peptide);
				
				if (entry==null) {
					return;
				}

				if (fittingParams.isEliminatedPeptide(peptide)) {
					return;
				}
				
				TFloatArrayList actual=new TFloatArrayList();
				for (ScoredObject<String> scoredObject : expectedConcentrations) {
					String column=scoredObject.y;
					String value = row.get(column);
					float concentration=Float.parseFloat(value);
					actual.add(concentration);
				}
				
				float[] actualArray = actual.toArray();
				
				if (!fittingParams.isTargetedProtein(protein)) {
					float mean = General.mean(actualArray);
					float cv=General.stdev(actualArray)/mean;
					for (int i = 0; i < bestIntensities.length; i++) {
						float rtInSec = rtAlignment.getAlignedRTInSec(entry);
						if (subRanges.get(i).contains(rtInSec)) {
							if (cv<fittingParams.getMinCVForAnchors()) {
								if (mean>bestIntensities[i]) {
									bestIntensities[i]=mean;
									bestAnchorPeptideModSeqs[i]=peptide;
								}
							} else if (bestIntensities[i]==0.0f&&cv<fittingParams.getMinCVForBadAnchors()) {
								if (mean>bestIntensitiesWithBadCVs[i]&&General.min(actualArray)>0.0f) {
									bestIntensitiesWithBadCVs[i]=mean;
									bestAnchorPeptideModSeqs[i]=peptide;
								}
							}
							break;
						}
					}
					// skip this for curve fitting
					return;
				}
			}
			
			public void cleanup() {
			}
		};
		
		if (dataFile.getName().toLowerCase().endsWith(".csv")) {
			TableParser.parseCSV(dataFile, muscle);
		} else {
			TableParser.parseTSV(dataFile, muscle);
		}

		return data;
	}

	protected static ArrayList<FitPeptide> fitCurves(final File outputDirectory, File dataFile, final ArrayList<ScoredObject<String>> expectedConcentrations,
			final float[] expected, AbstractDilutionCurveFittingParameters requiredAccessionText, boolean useLineNoise) throws FileNotFoundException, UnsupportedEncodingException {
		final PrintWriter reportWriter=new PrintWriter(new File(outputDirectory, "report.csv"), "UTF-8");
		reportWriter.println("peptide,protein,lod,loq,r2,m,b,noiseStdev,linearStdev,max");

		final ArrayList<FitPeptide> fitPeptides=new ArrayList<FitPeptide>();
		TableParserMuscle muscle = new TableParserMuscle() {
			public void processRow(Map<String, String> row) {
				String peptide = getPeptide(row);
				String protein = getProtein(row);
				
				if (requiredAccessionText!=null&&!requiredAccessionText.isTargetedProtein(protein)) {
					return;
				}

				if (requiredAccessionText.isEliminatedPeptide(peptide)) {
					return;
				}
				
				TFloatArrayList actual=new TFloatArrayList();
				for (ScoredObject<String> scoredObject : expectedConcentrations) {
					String column=scoredObject.y;
					try {
						float concentration=Float.parseFloat(row.get(column));
						actual.add(concentration);
					} catch (NullPointerException npe) {
						Logger.errorLine("Error parsing column ["+column+"], no column found!");
					} catch (NumberFormatException nfe) {
						Logger.errorLine("Error parsing column ["+column+"], value ["+row.get(column)+"] is not a number!");
					}
				}
				
				float[] actualArray = actual.toArray();

				float maxMeasuredValue = General.max(actualArray);
				actualArray=General.divide(actualArray, maxMeasuredValue);
				Pair<DilutionFit, Float> pair=process(peptide, protein, expected, actualArray, maxMeasuredValue, useLineNoise);
				DilutionFit bestFit=pair.x;
				if (bestFit==null) {
					Logger.logLine("Failed to fit "+peptide);
					return;
				}

				float lod=bestFit.getLOD();
				float loq=bestFit.getLOQ();
				
				reportWriter.println(peptide+","+protein+","+lod+","+loq+","+pair.y+","+bestFit.m+","+bestFit.b+","+bestFit.noiseStdev+","+bestFit.linearStdev+","+bestFit.maxValue);
				if (Float.isFinite(loq)&&loq<0) {
					fitPeptides.add(new FitPeptide(peptide, protein, bestFit, expected, actualArray));
				}
			}
			
			public void cleanup() {
			}
		};
		
		if (dataFile.getName().toLowerCase().endsWith(".csv")) {
			TableParser.parseCSV(dataFile, muscle);
		} else {
			TableParser.parseTSV(dataFile, muscle);
		}
		
		reportWriter.flush();
		reportWriter.close();
		
		Collections.sort(fitPeptides);

		Logger.logLine("Fit "+fitPeptides.size()+" total peptides.");
		
		return fitPeptides;
	}

	protected static HashMap<String, Map<String, TObjectFloatHashMap<String>>> extractUnknowns(File dataFile, final Map<String, TObjectFloatHashMap<String>> unknownSamples, final String requiredAccessionText) throws FileNotFoundException, UnsupportedEncodingException {
		final HashMap<String, Map<String, TObjectFloatHashMap<String>>> unknowns=new HashMap<>();
		if (unknownSamples==null||unknownSamples.size()==0) return unknowns;
		
		TableParserMuscle muscle = new TableParserMuscle() {
			public void processRow(Map<String, String> row) {
				String peptide = getPeptide(row);
				String protein = getProtein(row);
				
				if (requiredAccessionText!=null&&protein.indexOf(requiredAccessionText)==-1) {
					return;
				}
				
				TreeMap<String, TObjectFloatHashMap<String>> thisPeptideUnknowns=new TreeMap<>();
				unknowns.put(peptide, thisPeptideUnknowns);
				
				for (Entry<String, TObjectFloatHashMap<String>> entry : unknownSamples.entrySet()) {
					String group=entry.getKey();
					TObjectFloatHashMap<String> groupUnknowns=new TObjectFloatHashMap<>();
					thisPeptideUnknowns.put(group, groupUnknowns);

					entry.getValue().forEachEntry(new TObjectFloatProcedure<String>() {
						@Override
						public boolean execute(String column, float normalization) {
							float concentration;
							try {
								concentration=Float.parseFloat(row.get(column));
							} catch (NullPointerException npe) {
								Logger.errorLine("Failure to parse number from ["+row.get(column)+"] for column ["+column+"]");
								concentration=0.0f;
							} catch (NumberFormatException nfe) {
								Logger.errorLine("Failure to parse number from ["+row.get(column)+"] for column ["+column+"]");
								concentration=0.0f;
							}
							groupUnknowns.put(column, concentration*normalization);
							return true;
						}
					});
				}
			}
			
			public void cleanup() {
			}
		};
		
		if (dataFile.getName().toLowerCase().endsWith(".csv")) {
			TableParser.parseCSV(dataFile, muscle);
		} else {
			TableParser.parseTSV(dataFile, muscle);
		}

		Logger.logLine("Found unknowns for "+unknowns.size()+" total peptides.");
		
		return unknowns;
	}

	private static ArrayList<LibraryEntry> scheduleAssay(final File outputDirectory, final File targetDirectory,
			final File nontargetDirectory, AbstractDilutionCurveFittingParameters fittingParams,
			final HashMap<String, LibraryEntry> libraryEntryByPeptideModSeq, final AlignmentWithAnchors rtAlignment,
			Range rtInSecRange, final ArrayList<Range> subRanges, final String[] bestAnchorPeptideModSeqs,
			final float[] bestIntensities, final ArrayList<FitPeptide> fitPeptides)
			throws FileNotFoundException, UnsupportedEncodingException {
		ArrayList<LibraryEntry> targetEntries=new ArrayList<LibraryEntry>();
		boolean hitMaxDensity=false;
		float[] assayRT=new float[Math.round(rtInSecRange.getStop()+fittingParams.getWindowInMin()*60f)]; // N+W minutes in second increments
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

				assayDensity=incrementDensity(rtInSec, fittingParams.getWindowInMin(), assayDensity);
				addPeptideToAssay(assayWriter, entry, rtInSec, fittingParams.getWindowInMin());
				Logger.logLine("Using "+entry.getPeptideModSeq()+" from "+PSMData.accessionsToString(entry.getAccessions())+" as anchor (rt: "+(rtInSec/60f)+" mins, intensity: "+bestIntensities[i]+" for the RT range from "+(subRanges.get(i).getStart()/60f)+" min to "+(subRanges.get(i).getStop()/60f)+" min");
			} else {
				Logger.logLine("Failed to find good anchor for the RT range from "+(subRanges.get(i).getStart()/60f)+" min to "+(subRanges.get(i).getStop()/60f)+" min");
			}
		}
		
		int count=0;
		HashMap<String, ArrayList<FitPeptide>> targetPeptidesByProtein=new HashMap<String, ArrayList<FitPeptide>>();
		ArrayList<FitPeptide> nontargetedPeptides=new ArrayList<FitPeptide>();
		// assumes fitPeptides are sorted
		addpeptides:for (FitPeptide fit : fitPeptides) {
			ArrayList<FitPeptide> list=targetPeptidesByProtein.get(fit.proteinKey);
			if (list==null) {
				list=new ArrayList<DilutionCurveFitter.FitPeptide>();
				targetPeptidesByProtein.put(fit.proteinKey, list);
			}

			boolean keep=true;
			if (count<fittingParams.getTargetTotalNumberOfPeptides()) {
				if (list.size()<fittingParams.getMaxNumberPeptidesPerProtein()) {
					LibraryEntry entry=libraryEntryByPeptideModSeq.get(fit.peptideModSeq);
					if (entry==null) {
						System.out.println("NULL: "+fit.peptideModSeq);
						continue;
					}

					float rtInSec = rtAlignment.getAlignedRTInSec(entry);
					float[] testDensity=incrementDensity(rtInSec, fittingParams.getWindowInMin(), assayDensity);
					for (int i = 0; i < testDensity.length; i++) {
						if (testDensity[i]>fittingParams.getAssayMaxDensity()) {
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
				addPeptideToAssay(assayWriter, entry, rtInSec, fittingParams.getWindowInMin());
				targetEntries.add(entry.updateRetentionTime(rtInSec));

				ChartPanel panel=graph(fit.peptideModSeq, fit.expectedRelativeIntensities, fit.actualRelativeIntensities, fit.bestFit, Optional.empty());
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
			ChartPanel panel=graph(fit.peptideModSeq, fit.expectedRelativeIntensities, fit.actualRelativeIntensities, fit.bestFit, Optional.empty());
			Charter.writeAsPDF(panel.getChart(), new File(nontargetDirectory, fit.peptideModSeq+".pdf"), new Dimension(300, 300));
		}
		
		assayWriter.flush();
		assayWriter.close();
		Logger.logLine("Finished writing assay for "+targetPeptidesByProtein.size()+" proteins using "+count+" total peptides ("+numSingletons+" single peptide targets)");
		writeSchedulingGraph(outputDirectory, assayRT, assayDensity);
		return targetEntries;
	}

	static void writeLibraryEntries(SearchParameters params, final File exportLibraryFile,
			ArrayList<LibraryEntry> targetEntries) throws IOException, SQLException {
		LibraryFile exportLibrary=new LibraryFile();
		exportLibrary.openFile();
		exportLibrary.dropIndices();
		exportLibrary.addEntries(targetEntries);
		exportLibrary.addProteinsFromEntries(targetEntries);
		exportLibrary.addMetadata(params.toParameterMap());
		exportLibrary.createIndices();
		exportLibrary.saveAsFile(exportLibraryFile);
		exportLibrary.close();
	}

	static void writeSchedulingGraph(final File outputDirectory, float[] assayRT, float[] assayDensity) {
		XYTrace trace=new XYTrace(assayRT, assayDensity, GraphType.area, "Scheduling density");
		ChartPanel panel=Charter.getChart("Retention Time (min)", "Number of Peptides", true, trace);
		Charter.writeAsPDF(panel.getChart(), new File(outputDirectory, "assay_density.pdf"), new Dimension(600, 300));
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
    
	static HashMap<String, LibraryEntry> getLibraryData(SearchParameters params, File rtAlignFile) throws IOException, SQLException, DataFormatException {
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

	protected static float[] adjustForZeroConcentrations(final ArrayList<ScoredObject<String>> expectedConcentrations) {
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

	protected static Pair<ArrayList<ScoredObject<String>>, Map<String, TObjectFloatHashMap<String>>> getExpectedConcentrationsFromCSV(File sampleOrganizationFile) {
		final ArrayList<ScoredObject<String>> expectedConcentrations=new ArrayList<ScoredObject<String>>();
		final TreeMap<String, TObjectFloatHashMap<String>> unknowns=new TreeMap<>();
		
		System.out.println("Reading "+sampleOrganizationFile.getName()+"...");
		TableParserMuscle muscle = new TableParserMuscle() {
			public void processRow(Map<String, String> row) {
				String name=row.get("filename");
				String concentrationString = row.get("concentration");
				String normalizationString = row.get("normalization");
				try {
					float concentration=Float.parseFloat(concentrationString);
					expectedConcentrations.add(new ScoredObject<String>(concentration, name));
				} catch (NumberFormatException nfe) {
					// not a number, so parse as an unknown group
					TObjectFloatHashMap<String> list=unknowns.get(concentrationString);
					if (list==null) {
						list=new TObjectFloatHashMap<>();
						unknowns.put(concentrationString, list);
					}
					
					float normalization=1.0f;
					if (normalizationString!=null&&normalizationString.length()>0) {
						try {
							normalization=Float.parseFloat(normalizationString);
						} catch (NumberFormatException nfe2) {
							Logger.errorLine("Failed to parse normalization constant from "+name+", defaulting to 1.0");
						}
					}
					list.put(name, normalization);
				}
			}
			
			public void cleanup() {
			}
		};
		
		if (sampleOrganizationFile.getName().toLowerCase().endsWith(".csv")) {
			TableParser.parseCSV(sampleOrganizationFile, muscle);
		} else {
			TableParser.parseTSV(sampleOrganizationFile, muscle);
		}
		return new Pair<>(expectedConcentrations, unknowns);
	}
	
	public static Pair<DilutionFit, Float> process(String peptide, String protein, float[] expected, float[] actual, float maxMeasuredValue, boolean useLineNoise) {
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
				// noise only contains non-zero values
				noise.add(loggedActual.get(j));
			}
			float[] noiseArray = noise.toArray();
			float noiseMean = General.mean(noiseArray);
			float noiseMax=General.max(noiseArray);
			
			// fit everything after the crossover to a line
			TFloatArrayList linearX=new TFloatArrayList();
			TFloatArrayList linearY=new TFloatArrayList();
			boolean valuesAboveNoise=false;
			float sumVariance=0;
			int n=0;
			for (int j = crossOver+1; j < loggedActual.size(); j++) {
				n++;
				linearX.add(loggedExpected.get(j));
				linearY.add(loggedActual.get(j));
				float delta=loggedExpected.get(j)-loggedActual.get(j);
				sumVariance+=delta*delta;
				if (loggedActual.get(j)<noiseMean) {
					valuesAboveNoise=true;
				}
			}
			if (valuesAboveNoise||n==0) continue;
			
			float linearStdev=useLineNoise?(float)Math.sqrt(sumVariance/n):0.0f;
			float noiseStdev = General.stdev(noiseArray);

			// calculate equations
			Pair<Float, Float> equation=LinearRegression.getRegression(linearX.toArray(), linearY.toArray());
			
			// calculate deviation to find the best fit
			float rsquared=0;
			for (int j = 0; j < loggedExpected.size(); j++) {
				float x=loggedExpected.get(j);
				float actualY=loggedActual.get(j);

				float predictedY=equation.x*x+equation.y;
				if (predictedY<noiseMax) {
					predictedY=noiseMax;
				}
				
				float residual=actualY-predictedY;
				rsquared+=residual*residual;
			}
			
			fit=new DilutionFit(noiseMean, noiseMax, noiseStdev, linearStdev, equation.x, equation.y, lastZero, firstNonZero, maxMeasuredValue, rsquared);
			
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
				fit=new DilutionFit(noiseMean, noiseMax, noiseStdev, linearStdev, equation.x, equation.y, lastZero, firstNonZero, maxMeasuredValue, rsquared);
			}
//			if (peptide.equals("AGVLGALALGR")) {
//				System.out.println(crossOver+") "+rsquared+" ("+noise.size()+"/"+linearX.size()+") --> m:"+equation.x+", b:"+equation.y+", lastZero:"+lastZero); // FIXME
//			}
			
			// slope has to be at least 0.5
			if (fit.m>=0.5f&&rsquared<minRSquared) {
				minRSquared=rsquared;
				bestFit=fit;
			}
		}
		if (bestFit==null) bestFit=fit; // no good matches, so use last fit
		return new Pair<DilutionFit, Float>(bestFit,minRSquared);
	}
	
	public static ChartPanel graph(String peptide, float[] expected, float[] actual, DilutionFit bestFit, Optional<Map<String, TObjectFloatHashMap<String>>> unknowns) {
		expected=expected.clone();
		actual=actual.clone();
		
		float minNonZeroExpected=Float.MAX_VALUE;
		float minNonZeroActual=Float.MAX_VALUE;
		TFloatArrayList expectedFound=new TFloatArrayList();
		TFloatArrayList expectedMissing=new TFloatArrayList();
		TFloatArrayList actualFound=new TFloatArrayList();
		TFloatArrayList actualMissing=new TFloatArrayList();
		for (int i = 0; i < actual.length; i++) {
			if (actual[i]>0&&expected[i]>0) {
				if (actual[i]<minNonZeroActual) minNonZeroActual=actual[i];
				if (expected[i]<minNonZeroExpected) minNonZeroExpected=expected[i];
			}
		}
		
		for (int i = 0; i < actual.length; i++) {
			if (actual[i]<=0&&expected[i]<=0) {
				actualMissing.add(minNonZeroActual/10f);
				expectedMissing.add(minNonZeroExpected/10f);
			} else if (actual[i]<=0) {
				actualMissing.add(minNonZeroActual/10f);
				expectedMissing.add(expected[i]);
			} else if (expected[i]<=0) {
				actualMissing.add(actual[i]);
				expectedMissing.add(minNonZeroExpected/10f);
			} else {
				actualFound.add(actual[i]);
				expectedFound.add(expected[i]);
			}
			if (actual[i]<=0) actual[i]=minNonZeroActual/10f;
			if (expected[i]<=0) expected[i]=minNonZeroExpected/10f;

		}
		
		float lod=(float)Math.pow(10, bestFit.getLOD());
		float loq=(float)Math.pow(10, bestFit.getLOQ());
		float maxExpectedWithMargin=General.max(expected)*10;
		float minExpected = General.min(expected);
		float minExpectedWithMargin=minExpected/10;

		float adjustmentForUnknowns=General.max(expected);
		if (unknowns.isPresent()) {
			adjustmentForUnknowns=(float)Math.pow(10, unknowns.get().size());
		}
		
		XYTrace lodTrace=new XYTrace(new float[] {lod}, new float[] {bestFit.getUnloggedPredicted(lod)}, GraphType.bighollowpoint, "LOD="+lod, Color.gray, 10f);
		XYTrace loqTrace=new XYTrace(new float[] {minExpected/5, maxExpectedWithMargin*adjustmentForUnknowns, Float.NaN, loq, loq}, new float[] {loq, loq, Float.NaN, minExpectedWithMargin, maxExpectedWithMargin/10}, GraphType.dashedline, "LOQ="+loq, Color.red, 2f);
		
		XYTrace actualTrace=new XYTrace(expectedFound.toArray(), actualFound.toArray(), GraphType.bigpoint, peptide, Color.BLACK, 10f);
		XYTrace actualMissingTrace=new XYTrace(expectedMissing.toArray(), actualMissing.toArray(), GraphType.bighollowpoint, "Missing", Color.BLACK, 10f);

		ChartPanel panel=Charter.getChart("Expected", "Actual", true, actualTrace, actualMissingTrace, lodTrace);

		XYPlot plot = panel.getChart().getXYPlot();
		int currentCount=plot.getDatasetCount();

		// DRAW BOX PLOT FOR UNKNOWNS
		int boxItemWidth = 12;
		int unknownCount=0;
		TreeMap<Double, String> unknownXLocations=new TreeMap<>();
		if (unknowns.isPresent()) {
			float currentXValue=maxExpectedWithMargin/3;
			Map<String, TObjectFloatHashMap<String>> unknownMap=unknowns.get();
			for (Entry<String, TObjectFloatHashMap<String>> entry : unknownMap.entrySet()) {
				TObjectFloatHashMap<String> unknownData=entry.getValue();
				currentXValue=currentXValue*4.0f;
				unknownXLocations.put(Double.valueOf(currentXValue), entry.getKey());
				
				NumberBoxAndWhiskerXYDataset boxplotDataset=new NumberBoxAndWhiskerXYDataset("Unknowns");
				float[] values = General.divide(unknownData.values(), bestFit.maxValue);
				BoxAndWhiskerItem stats = Boxplotter.calculateINFProtectedBoxAndWhiskerStatistics(values, maxExpectedWithMargin, minExpectedWithMargin, minExpectedWithMargin, true);
				
				boxplotDataset.add(currentXValue, stats);
				
				plot.setDataset(currentCount, boxplotDataset);
				XYBoxPlotterRenderer boxplotRenderer = new XYBoxPlotterRenderer(boxItemWidth, true);
				boxplotRenderer.setSeriesVisibleInLegend(0, false);
				boxplotRenderer.setSeriesPaint(0, colors[unknownCount%colors.length]);
				plot.setRenderer(currentCount, boxplotRenderer);
				unknownCount++;
				currentCount++;
			}
		}
		
		// DRAW LOQ
		XYGraphingTrace loqGraphingTrace=new XYGraphingTrace(loqTrace);
		XYSeriesCollection loqDataset=new XYSeriesCollection();
		loqDataset.addSeries(loqGraphingTrace.getSeries());
		plot.setDataset(currentCount, loqDataset);
		plot.setRenderer(currentCount, loqGraphingTrace.getRenderer());
		currentCount++;
		
		// DRAW STANDARD DEVIATION SHADING
		TFloatArrayList expectedPlusLODList=new TFloatArrayList(expected);
		expectedPlusLODList.add(lod);
		expectedPlusLODList.sort();
		float[] expectedPlusLOD=expectedPlusLODList.toArray();
		float[] predicted = bestFit.getUnloggedPredicted(expectedPlusLOD);
		
		DefaultIntervalXYDataset dataset=new DefaultIntervalXYDataset();
		float[] expectedLower=bestFit.getUnloggedLowerError(expectedPlusLOD);
		float[] expectedUpper=bestFit.getUnloggedUpperError(expectedPlusLOD);
		dataset.addSeries("Calculated", General.toDoubleArray(new float[][] {expectedPlusLOD, expectedPlusLOD, expectedPlusLOD, predicted, expectedLower, expectedUpper}));
		plot.setDataset(currentCount, dataset);
		
		DeviationRenderer renderer = new DeviationRenderer(true, false);
		renderer.setSeriesFillPaint(0, new Color(255, 255, 0, 175));
		renderer.setSeriesStroke(0, new BasicStroke(2, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 0.0f, new float[] {3.0f, 5.0f}, 0.0f));
		renderer.setSeriesPaint(0, Color.gray);
		renderer.setSeriesVisibleInLegend(0, false);
		plot.setRenderer(currentCount, renderer);
		
		ValueAxis domain=plot.getDomainAxis();
		ValueAxis range=plot.getRangeAxis();
		ExtendedLogAxis newDomain = new ExtendedLogAxis("Expected") {
			@Override
			protected LogTick getMinorTick(TextAnchor textAnchor, double v) {
				if (v>1) return null;
				return super.getMinorTick(textAnchor, v);
			}
			@Override
			protected LogTick getMajorTick(TextAnchor textAnchor, double v) {
				if (v>1) return null;
				return super.getMajorTick(textAnchor, v);
			}
		};
		for (Entry<Double, String> entry : unknownXLocations.entrySet()) {
			newDomain.addAutomaticTicks(entry.getKey(), entry.getValue());
		}
		newDomain.setLabelFont(domain.getLabelFont());
		newDomain.setTickLabelFont(domain.getTickLabelFont());
		
		ExtendedLogAxis newRange =new ExtendedLogAxis("Actual") {
			@Override
			protected LogTick getMinorTick(TextAnchor textAnchor, double v) {
				if (v>1) return null;
				return super.getMinorTick(textAnchor, v);
			}
			@Override
			protected LogTick getMajorTick(TextAnchor textAnchor, double v) {
				if (v>1) return null;
				return super.getMajorTick(textAnchor, v);
			}
		};
		
		newRange.setLabelFont(range.getLabelFont());
		newRange.setTickLabelFont(range.getTickLabelFont());

		plot.setDomainAxis(newDomain);
		plot.setRangeAxis(newRange);
		return panel;
	}
	
	private static String getPeptide(Map<String, String> row) {
		String peptide=row.get("Peptide"); // from EncyclopeDIA
		if (peptide==null) peptide=row.get("Peptide Modified Sequence"); // from Skyline
		if (peptide==null) peptide=row.get("Peptide Modified Sequence Monoisotopic Masses"); // from Skyline
		return peptide;
	}

	private static String getProtein(Map<String, String> row) {
		String protein=row.get("Protein"); // from EncyclopeDIA
		if (protein==null) protein=row.get("Protein Name"); // from Skyline
		if (protein==null) protein="Unknown Protein";
		return protein;
	}

	public static Color[] colors=new Color[] {new Color(95,158,160), new Color(240,230,140), new Color(255, 215, 0), new Color(205,173,0), new Color(139,117,0), new Color(139,76,57), new Color(139,0,0)};

	public static class AlignmentWithAnchors {
		final RetentionTimeFilter rtAlignment;
		final TObjectFloatHashMap<String> knownRTInSecs;
		public AlignmentWithAnchors(RetentionTimeFilter rtAlignment, TObjectFloatHashMap<String> knownRTInSecs) {
			this.rtAlignment = rtAlignment;
			this.knownRTInSecs = knownRTInSecs;
		}

		public float getAlignedRTInSec(LibraryEntry entry) {
			return getAlignedRTInSec(entry, true);
		}
		public boolean isKnown(LibraryEntry entry) {
			return knownRTInSecs.contains(entry.getPeptideModSeq());
		}
		
		public float getAlignedRTInSec(LibraryEntry entry, boolean quiet) {
			if (knownRTInSecs.contains(entry.getPeptideModSeq())) {
				return knownRTInSecs.get(entry.getPeptideModSeq());
			}
			if (!quiet)	Logger.errorLine("Potential problem: had to look up "+entry.getPeptideModSeq());
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

		public float[] getExpectedRelativeIntensities() {
			return expectedRelativeIntensities;
		}

		public void setExpectedRelativeIntensities(float[] expectedRelativeIntensities) {
			this.expectedRelativeIntensities = expectedRelativeIntensities;
		}

		public float[] getActualRelativeIntensities() {
			return actualRelativeIntensities;
		}

		public void setActualRelativeIntensities(float[] actualRelativeIntensities) {
			this.actualRelativeIntensities = actualRelativeIntensities;
		}

		public String getPeptideModSeq() {
			return peptideModSeq;
		}

		public String getProteinKey() {
			return proteinKey;
		}

		public DilutionFit getBestFit() {
			return bestFit;
		}
		
	}
	
	private static final float NUM_STDEVS_FOR_LOQ=3.0f;
	
	public static class DilutionFit {
		private final float noiseMean;
		private final float noiseMax;
		private final float noiseStdev;
		private final float linearStdev;
		private final float m;
		private final float b;
		private final float lastZero;
		private final float firstNonZero;
		private final float maxValue;
		private final float r2;
		
		public DilutionFit(float noiseMean, float noiseMax, float noiseStdev, float linearStdev, float m, float b, float lastZero, float firstNonZero, float maxValue, float r2) {
			this.noiseMean = noiseMean;
			this.noiseMax=noiseMax;
			this.noiseStdev = noiseStdev;
			this.linearStdev=linearStdev;
			this.m = m;
			this.b = b;
			this.lastZero=lastZero;
			this.firstNonZero=firstNonZero;
			this.maxValue=maxValue;
			this.r2=r2;
		}
		
		public float getPredicted(float x) {
			float expectedY=m*x+b;
			if (expectedY<noiseMax) {
				return noiseMax;
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
		public float[] getUnloggedUpperError(float[] xs) {
			float loq=getLOQ();
			float[] expectedYs=new float[xs.length];
			for (int i = 0; i < expectedYs.length; i++) {
				float x=xs[i];

				float loggedX=Log.log10(x);
				float expectedY=getPredicted(loggedX)+NUM_STDEVS_FOR_LOQ*linearStdev;
				if (expectedY<loq) {
					expectedY=loq;
				}
				expectedYs[i]=(float)Math.pow(10, expectedY);
			}
			return expectedYs;
		}

		public float[] getUnloggedLowerError(float[] xs) {
			float[] expectedYs=new float[xs.length];
			for (int i = 0; i < expectedYs.length; i++) {
				float x=xs[i];

				float loggedX=Log.log10(x);
				float expectedY=getPredicted(loggedX);
				if (expectedY>noiseMax) {
					expectedY-=NUM_STDEVS_FOR_LOQ*linearStdev;
				} else {
					// noise level
					expectedY-=NUM_STDEVS_FOR_LOQ*getStdev();
				}
				expectedYs[i]=(float)Math.pow(10, expectedY);
			}
			return expectedYs;
		}
		
		public float getLOD() {
			//noiseValue=mx+b
			if (m==0) return 0;
			if (r2>1000) return 0;
			float lod = Math.max(lastZero, Math.min(0f, (noiseMax-b)/m));
			if (Float.isInfinite(lod)) return 0;
			return lod;
		}
		public float getLOQ() {
			if (m==0) return 0;
			if (r2>1000) return 0;
			float target=noiseMax+NUM_STDEVS_FOR_LOQ*getStdev();
			float loq = Math.max(firstNonZero, Math.min(0f, (target-b)/m));
			if (Float.isInfinite(loq)) return 0;
			return loq;
		}
		
		public float getStdev() {
			return Math.max(linearStdev, noiseStdev);
		}
		
		public float getMaxValue() {
			return maxValue;
		}
	}
}
