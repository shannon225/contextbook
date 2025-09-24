package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.LibraryReportExtractor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.*;
import edu.washington.gs.maccoss.encyclopedia.filereaders.*;
import edu.washington.gs.maccoss.encyclopedia.gui.general.*;
import edu.washington.gs.maccoss.encyclopedia.gui.massspec.ChromatogramCharter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.StringUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.*;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.CombinedRangeXYPlot;
import org.jfree.ui.RectangleInsets;
import org.jfree.ui.TextAnchor;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.text.AttributedString;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.zip.DataFormatException;

public class MultiResultsBrowserPanel extends JPanel {
	private static final long serialVersionUID=1L;

	private static final int RT_EXTRACTION_MARGIN_IN_SEC=15;
	
	private final FileChooserPanel elibFileChooser;
	private final FileChooserPanel libraryFileChooser;
	private final JSplitPane split=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	private final JSplitPane split1=new JSplitPane(JSplitPane.VERTICAL_SPLIT);
	private final JSplitPane split2=new JSplitPane(JSplitPane.VERTICAL_SPLIT);
	
	private final JTable sampleTable;
	private final TableRowSorter<TableModel> sampleRowSorter;
	private final SampleTableModel sampleModel;
	private final JTable peptideTable;
	private final TableRowSorter<TableModel> peptideRowSorter;
	private final JTextField jtfFilter;
	private final JCheckBox jtfNotFilter=new JCheckBox("NOT");
	private final MultiPeptideResultsTableModel peptideModel;
	private final SearchParameters parameters;
	private final ChartPanel barChart;
	private final ChartPanel stackedBarChart;
	private final JComboBox<Integer> minimumNumberOfFragments=new JComboBox<Integer>(new Integer[] {0, 1, 2, 3, 4, 5});
	private final JComboBox<Integer> numberOfColumns=new JComboBox<Integer>(new Integer[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
	private final JCheckBox simplifyPlots=new JCheckBox("Simplify plots");
	private final JCheckBox precursorsPlots=new JCheckBox("Precursors");
	private final JCheckBox fragmentsPlots=new JCheckBox("Fragments");
	private final JCheckBox sgSmoothBox=new JCheckBox("Smooth");
	private final JCheckBox backgroundSubtractBox=new JCheckBox("Background Subtract");
	
	private final int defaultMinimumNumberOfFragmentsIndex=3;
	private final int defaultNumberOfColumnsIndex=1;
	private volatile LibraryFile library; // CAN BE NULL, CAN BE MUTATED!
	
	public MultiResultsBrowserPanel(SearchParameters parameters) {
		super(new BorderLayout());
		
		this.parameters=parameters;

		sampleModel=new SampleTableModel();
		sampleTable=new JTable(sampleModel) {
			private static final long serialVersionUID=1L;

			@Override
			public Object getValueAt(int row, int column) {
				if (column==0) return row+1;
				return super.getValueAt(row, column);
			}
		};
		sampleRowSorter=new TableRowSorter<TableModel>(sampleTable.getModel());
		sampleTable.setRowSorter(sampleRowSorter);

		peptideModel=new MultiPeptideResultsTableModel();
		peptideTable=new JTable(peptideModel) {

			@Override
			public Object getValueAt(int row, int column) {
				if (column==0) return row+1;
				return super.getValueAt(row, column);
			}
		};
		peptideRowSorter=new TableRowSorter<TableModel>(peptideTable.getModel());
		peptideTable.setRowSorter(peptideRowSorter);

		jtfFilter=new JTextField();
		jtfFilter.getDocument().addDocumentListener(new DocumentListener() {
			@Override
			public void insertUpdate(DocumentEvent e) {
				updateFilter();
			}

			@Override
			public void removeUpdate(DocumentEvent e) {
				updateFilter();
			}

			@Override
			public void changedUpdate(DocumentEvent e) {
				throw new UnsupportedOperationException("Not supported yet.");
			}
		});
		jtfNotFilter.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				updateFilter();
			}
		});

		peptideTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			@Override
			public void valueChanged(ListSelectionEvent e) {
				updateToSelectedPeptide();
			}
		});

		libraryFileChooser=new FileChooserPanel(null, "Library", new SimpleFilenameFilter(LibraryFile.DLIB, LibraryFile.ELIB), false) {
			private static final long serialVersionUID=1L;

			@Override
			public void update(File... filenames) {
				super.update(filenames);
				if (filenames!=null&&filenames.length>0&&filenames[0]!=null) {
					updateLibrary(filenames[0]);
				}
			}
		};

		elibFileChooser=new FileChooserPanel(null, "Results", new SimpleFilenameFilter(LibraryFile.DLIB, LibraryFile.ELIB), true) {
			private static final long serialVersionUID=1L;

			@Override
			public void update(File... filenames) {
				super.update(filenames);
				if (filenames!=null&&filenames.length>0&&filenames[0]!=null) {
					updateTables(filenames[0]);
				}
			}
		};
		numberOfColumns.setSelectedIndex(defaultNumberOfColumnsIndex);
		numberOfColumns.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				updateToSelectedPeptide();
			}
		});
		
		minimumNumberOfFragments.setSelectedIndex(defaultMinimumNumberOfFragmentsIndex);
		Integer minimumNumberOfTransitions=(Integer)minimumNumberOfFragments.getSelectedItem();
		peptideModel.filterTable(minimumNumberOfTransitions);
		
		minimumNumberOfFragments.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				Integer minimumNumberOfTransitions=(Integer)minimumNumberOfFragments.getSelectedItem();
				peptideModel.filterTable(minimumNumberOfTransitions);
			}
		});

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(elibFileChooser);
		options.add(libraryFileChooser);
		options.add(new LabeledComponent("Minimum # Fragments", minimumNumberOfFragments));
		options.add(new LabeledComponent("Number of Columns", numberOfColumns));
		

		JPanel checkBoxes=new JPanel();
		checkBoxes.setLayout(new BoxLayout(checkBoxes, BoxLayout.LINE_AXIS));
		simplifyPlots.setBackground(LabeledComponent.BACKGROUND_COLOR);
		checkBoxes.add(simplifyPlots);
		options.add(checkBoxes);
		simplifyPlots.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				updateToSelectedPeptide();
			}
		});
		precursorsPlots.setBackground(LabeledComponent.BACKGROUND_COLOR);
		precursorsPlots.setSelected(true);
		checkBoxes.add(precursorsPlots);
		options.add(checkBoxes);
		precursorsPlots.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				updateToSelectedPeptide();
			}
		});
		fragmentsPlots.setBackground(LabeledComponent.BACKGROUND_COLOR);
		fragmentsPlots.setSelected(true);
		checkBoxes.add(fragmentsPlots);
		
		sgSmoothBox.setSelected(true);
		sgSmoothBox.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				updateToSelectedPeptide();
			}
		});

		backgroundSubtractBox.setSelected(true);
		backgroundSubtractBox.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				updateToSelectedPeptide();
			}
		});
		checkBoxes.add(sgSmoothBox);
		checkBoxes.add(backgroundSubtractBox);
		
		options.add(checkBoxes);
		fragmentsPlots.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				updateToSelectedPeptide();
			}
		});
		
		split2.setBorder(null); // remove border of nested split so they seem in the same level
		split1.setTopComponent(new JScrollPane(sampleTable, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED));
		split1.setBottomComponent(split2);
		split2.setTopComponent(new JScrollPane(peptideTable, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED));
		
		barChart=getBarChart(new String[] {}, new float[] {});
		stackedBarChart=getBarChart(new String[] {}, new float[] {});

		JTabbedPane barChartTabs=new JTabbedPane();
		barChartTabs.addTab("Relative Intensities", stackedBarChart);
		barChartTabs.addTab("Total Intensities", barChart);
		split2.setBottomComponent(barChartTabs);

		JPanel searchPanel=new JPanel(new BorderLayout());
		searchPanel.add(new JLabel("Search:"), BorderLayout.WEST);
		searchPanel.add(jtfFilter, BorderLayout.CENTER);
		searchPanel.add(jtfNotFilter, BorderLayout.EAST);

		JPanel left=new JPanel(new BorderLayout());
		left.add(options, BorderLayout.NORTH);
		left.add(split1, BorderLayout.CENTER);
		left.add(searchPanel, BorderLayout.SOUTH);
		
		
		
		split.setLeftComponent(left);
		split.setRightComponent(new JLabel("Select a peptide!"));
		

		setLayout(new BorderLayout());
		add(split, BorderLayout.CENTER);
	}

	private void updateFilter() {
		String text=jtfFilter.getText();

		if (text.trim().length()==0) {
			peptideRowSorter.setRowFilter(null);
		} else if (jtfNotFilter.isSelected()) {
			peptideRowSorter.setRowFilter(RowFilter.notFilter(RowFilter.regexFilter("(?i)"+text)));
		} else {
			peptideRowSorter.setRowFilter(RowFilter.regexFilter("(?i)"+text));
		}
	}

	public ChartPanel getBarChart(String[] categories, float[] intensities) {
		return Charter.getBarChart(null, "Sample", "Total Intensity", categories, intensities);
	}
	
	public ChartPanel getStackedBarChart(PeptidePrecursor peptide, String[] categories, Spectrum[] dataArray) {
		Spectrum libraryEntry=null;
		if (library!=null) {
			try {
				ArrayList<LibraryEntry> entries=library.getEntries(peptide.getPeptideModSeq(), peptide.getPrecursorCharge(), false);
				if (entries.size()>0) {
					libraryEntry=entries.get(0);
				}
			} catch (Exception e) {
				Logger.errorLine("Error reading library entry!");
				Logger.errorException(e);
			}
		}
		
		return FragmentIonConsistencyCharter.getBarChart(peptide, Optional.ofNullable(libraryEntry), dataArray, categories, parameters);
	}

	public void askForResults() {
		elibFileChooser.askForFiles();
	}
	
	public void updateLibrary(final File f) {
		SwingWorkerProgress<LibraryFile> worker=new SwingWorkerProgress<LibraryFile>((Frame)SwingUtilities.getWindowAncestor(this), "Please wait...", "Reading Library") {
			@Override
			protected LibraryFile doInBackgroundForReal() throws Exception {
				LibraryFile.OPEN_IN_PLACE=true;
				LibraryInterface ilib=BlibToLibraryConverter.getFile(f);
				LibraryFile.OPEN_IN_PLACE=false;
				if (!(ilib instanceof LibraryFile)) {
					throw new EncyclopediaException("Sorry, can't load this type of library file "+ilib.getClass().getName());
				}
				LibraryFile library=(LibraryFile)ilib;
				return library;
			}
			@Override
			protected void doneForReal(LibraryFile t) {
				Logger.logLine("Finished loading library, updating GUI");
				library=t;
				updateToSelectedPeptide();
			}
		};
		worker.execute();
	}
	
	public void updateTables(final File f) {
		SwingWorkerProgress<Pair<ArrayList<StripeFileInterface>, ArrayList<PeptideReportData>>> worker=new SwingWorkerProgress<Pair<ArrayList<StripeFileInterface>, ArrayList<PeptideReportData>>>((Frame)SwingUtilities.getWindowAncestor(this), "Please wait...", "Reading Library") {
			@Override
			protected Pair<ArrayList<StripeFileInterface>, ArrayList<PeptideReportData>> doInBackgroundForReal() throws Exception {
				LibraryFile.OPEN_IN_PLACE=true;
				LibraryInterface ilib=BlibToLibraryConverter.getFile(f);
				LibraryFile.OPEN_IN_PLACE=false;
				if (!(ilib instanceof LibraryFile)) {
					throw new EncyclopediaException("Sorry, can't load this type of library file "+ilib.getClass().getName());
				}
				LibraryFile library=(LibraryFile)ilib;
				
				Pair<ArrayList<String>, ArrayList<PeptideReportData>> pair=LibraryReportExtractor.extractMatrix(library, parameters.getAAConstants());
				
				ArrayList<String> sampleNames=pair.x;
				ArrayList<StripeFileInterface> stripeFiles=new ArrayList<StripeFileInterface>();

				for (String sampleName : sampleNames) {
					Logger.logLine("Trying to load "+sampleName);
					StripeFileInterface file=StripeFileGenerator.getFile(new File(f.getParentFile(), sampleName), parameters, true);
					stripeFiles.add(file);
				}
				
				return new Pair<ArrayList<StripeFileInterface>, ArrayList<PeptideReportData>>(stripeFiles, pair.y);
			}
			@Override
			protected void doneForReal(Pair<ArrayList<StripeFileInterface>, ArrayList<PeptideReportData>> t) {
				Logger.logLine("Finished loading data, updating GUI ("+t.x.size()+" files, "+t.y.size()+" peptides)");
				peptideModel.updateEntries(t.y);
				sampleModel.updateEntries(t.x);

				if (split1.getDividerLocation()<=40||split2.getDividerLocation()<=40) {
					split2.setDividerLocation(200);
					split1.setDividerLocation(200);
				}
				
				// roughly 2 rows is a good heuristic for the GUI
				int numOfSamples=sampleModel.getRows().size();
				int numExpectedColumns=Math.round(numOfSamples/2.0f);
				if (numExpectedColumns>5) numExpectedColumns=5; // after 6 columns the plot gets messy
				
				// this updates the GUI plots 
				numberOfColumns.setSelectedIndex(numExpectedColumns);
			}
		};
		worker.execute();
	}
	
	public void updateToSelectedPeptide() {
		int[] selection=peptideTable.getSelectedRows();
		if (selection.length<=0) return;
		
		PeptideReportData entry=peptideModel.getSelectedRow(peptideTable.convertRowIndexToModel(selection[0]));
		resetPeptide(entry, sampleModel.getRows());
	}

	private static final DecimalFormat formatter=new DecimalFormat("0.#E0");
	public void resetPeptide(PeptideReportData entry, ArrayList<StripeFileInterface> files) {
		int location=split.getDividerLocation();
		if (location<=5) {
			location=200;
		}
		
		String[] origSampleNames=new String[files.size()];
		float[] totalTICs=new float[files.size()];
		Spectrum[] dataArray=new QuantitativeDIAData[files.size()];
		for (int i=0; i<origSampleNames.length; i++) {
			origSampleNames[i]=files.get(i).getOriginalFileName();
			dataArray[i]=entry.getQuantitativeData(origSampleNames[i]);
			if (dataArray[i]!=null) {
				totalTICs[i]=dataArray[i].getTIC();
			}
		}
		String[] sampleNames=StringUtils.getUniquePortion(origSampleNames);
		barChart.setChart(getBarChart(sampleNames, totalTICs).getChart());
		
		stackedBarChart.setChart(getStackedBarChart(entry, sampleNames, dataArray).getChart());
		
		int cols=(Integer)numberOfColumns.getSelectedItem();
		boolean isSimplify=simplifyPlots.isSelected();
		boolean isPrecursors=precursorsPlots.isSelected();
		boolean isFragments=fragmentsPlots.isSelected();
		if (files.size()<cols) cols=files.size();

		try {
			JPanel chromatogramsPanel=new JPanel(new GridLayout(0, isSimplify?1:cols));
			chromatogramsPanel.setBackground(Color.WHITE);
			FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
			FragmentIon[] primaryIonObjects=model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), false);
			Logger.logLine("Graphing "+entry.getPeptideModSeq()+" ("+primaryIonObjects.length+")"+"...");
			
			double globalMaxYFragment=0.0;
			double globalMaxYPrecursor=0.0;
			
			ArrayList<ArrayList<XYTrace>> allFragmentTraces=new ArrayList<ArrayList<XYTrace>>();
			ArrayList<ArrayList<XYTrace>> allPrecursorTraces=new ArrayList<ArrayList<XYTrace>>();
			for (int i=0; i<sampleNames.length; i++) {
				String sampleName = sampleNames[i];
				QuantitativeDIAData quantitativeData=entry.getQuantitativeData(sampleName);
				StripeFileInterface file=files.get(i);

				ArrayList<XYTrace> fragmentTraces;
				ArrayList<XYTrace> precursorTraces;
				if (quantitativeData==null) {
					fragmentTraces=new ArrayList<>();
					precursorTraces=new ArrayList<>();
				} else {
					precursorTraces=extractPrecursorTraces(sampleName, quantitativeData, file);
					fragmentTraces=extractFragmentTraces(primaryIonObjects, sampleName, quantitativeData, file);
				}

				if (fragmentTraces!=null) {
					for (XYTrace xyTrace : fragmentTraces) {
						if (xyTrace.getType()==GraphType.line&&xyTrace.size()>0) {
							globalMaxYFragment=Math.max(globalMaxYFragment, xyTrace.getMaxY());
						}
					}
				} else {
					Logger.logLine("Couldn't extract fragments from "+sampleName+" for "+entry.getPeptideModSeq());
				}
				allFragmentTraces.add(fragmentTraces);

				if (precursorTraces!=null) {
					for (XYTrace xyTrace : precursorTraces) {
						if (xyTrace.getType()==GraphType.line) {
							globalMaxYPrecursor=Math.max(globalMaxYPrecursor, xyTrace.getMaxY());
						}
					}
				} else {
					Logger.logLine("Couldn't extract precursors	 from "+sampleName+" for "+entry.getPeptideModSeq());
				}
				allPrecursorTraces.add(precursorTraces);
			}

			globalMaxYFragment=globalMaxYFragment*1.05;
			globalMaxYPrecursor=globalMaxYPrecursor*1.05;

			CombinedRangeXYPlot parent=null;
			for (int i=0; i<sampleNames.length; i++) {
				if (isSimplify) {
					ChartPanel fragmentChart=Charter.getChart("Retention Time (min)", "Intensity", false, globalMaxYFragment, allFragmentTraces.get(i).toArray(new XYTrace[0]));
					
					fragmentChart.getChart().getXYPlot().clearAnnotations();
					
					ValueAxis domainAxis = fragmentChart.getChart().getXYPlot().getDomainAxis();
					ValueAxis rangeAxis = fragmentChart.getChart().getXYPlot().getRangeAxis();
					XYTextAnnotation annotation = new XYTextAnnotation(sampleNames[i], domainAxis.getLowerBound(), rangeAxis.getUpperBound());
					annotation.setTextAnchor(TextAnchor.TOP_LEFT);
					annotation.setFont(new Font("News Gothic MT", Font.BOLD, 14));
					//fragmentChart.getChart().getXYPlot().addAnnotation(annotation);
					
					domainAxis.setLabel(null);
					domainAxis.setTickLabelFont(new Font("News Gothic MT", Font.PLAIN, 12));
					
					if (i%cols==0) {
						// ADD label the domain of the left most plots (FIXME FIND A BETTER WAY TO DO THIS)
				        parent = new CombinedRangeXYPlot(rangeAxis);

						JFreeChart chart = new JFreeChart(parent);
						chart.setPadding(new RectangleInsets(10, 10, 10, 10));
						final ChartPanel chartPanel=new ChartPanel(chart, false);
						chartPanel.getChart().removeLegend();
						chartPanel.getChart().setBackgroundPaint(Color.white);
						chartPanel.setMinimumDrawWidth(0);
						chartPanel.setMinimumDrawHeight(0);
						chartPanel.setMaximumDrawWidth(Integer.MAX_VALUE);
						chartPanel.setMaximumDrawHeight(Integer.MAX_VALUE);
				        chromatogramsPanel.add(chartPanel);
				        
						rangeAxis.setLabelFont(new Font("News Gothic MT", Font.PLAIN, 12));
					} else {
						rangeAxis.setAttributedLabel((AttributedString)null);
						rangeAxis.setLabel(null);
						rangeAxis.setTickLabelsVisible(false);
						rangeAxis.setTickMarksVisible(true);
						rangeAxis.setLabelFont(new Font("News Gothic MT", Font.PLAIN, 12));
						
					}
					
					if (parent!=null) {
						parent.add(fragmentChart.getChart().getXYPlot(), 1);
					}
				} else {
					// !simplify
					Optional<ArrayList<XYTrace>> precursors=isPrecursors?Optional.ofNullable(allPrecursorTraces.get(i)):Optional.empty();
					Optional<ArrayList<XYTrace>> fragments=isFragments?Optional.ofNullable(allFragmentTraces.get(i)):Optional.empty();
					ChartPanel panel=ChromatogramCharter.createChart(precursors, fragments, globalMaxYPrecursor, globalMaxYFragment);
					
					panel.getChart().setTitle(sampleNames[i]);
					chromatogramsPanel.add(panel);
					
				}
			}
	
			split.setRightComponent(chromatogramsPanel);

		} catch (DataFormatException sqle) {
			Logger.errorLine("Error reading raw files!");
			Logger.errorException(sqle);
		} catch (SQLException sqle) {
			Logger.errorLine("Error reading raw files!");
			Logger.errorException(sqle);
		} catch (IOException ioe) {
			Logger.errorLine("Error reading raw files!");
			Logger.errorException(ioe);
		}
		
		split.setDividerLocation(location);
	}
	
	private ArrayList<XYTrace> extractPrecursorTraces(String sampleName, QuantitativeDIAData quantitativeData, StripeFileInterface file) throws IOException, SQLException, DataFormatException {
		Range rangeInSec=quantitativeData.getRtScanRange();
		float minRTInSec = rangeInSec.getStart()-RT_EXTRACTION_MARGIN_IN_SEC;
		float maxRTInSec = rangeInSec.getStop()+RT_EXTRACTION_MARGIN_IN_SEC;
		ArrayList<PrecursorScan> precursors=file.getPrecursors(minRTInSec, maxRTInSec);
		ArrayList<PrecursorScan> trimmedPrecursors=new ArrayList<>();
		for (PrecursorScan spectrum : precursors) {
			if (quantitativeData.getPrecursorMZ()>spectrum.getIsolationWindowLower()&&quantitativeData.getPrecursorMZ()<spectrum.getIsolationWindowUpper()) {
				trimmedPrecursors.add(spectrum);
			}
		}
		precursors=trimmedPrecursors;
		XYTraceInterface[] traceArray=ChromatogramExtractor.extractPrecursorChromatograms(parameters.getPrecursorTolerance(), quantitativeData.getPrecursorMZ(), quantitativeData.getPrecursorCharge(), precursors, sgSmoothBox.isSelected(), backgroundSubtractBox.isSelected());

		double maxY=0.0;
		for (int i = 0; i < traceArray.length; i++) {
			maxY=Math.max(maxY, General.max(traceArray[i].toArrays().y));
		}
		
		ArrayList<XYTrace> traces=new ArrayList<>();
		for (int i = 0; i < traceArray.length; i++) {
			if (traceArray[i] instanceof XYTrace) {
				traces.add((XYTrace)traceArray[i]);
			}
		}
		
		if (traces.size()>0) {
			// extra 0 point in case there is no data shown (or all 0s)
			traces.add(new XYTrace(new double[] {rangeInSec.getStart()/60f-Float.MIN_VALUE, rangeInSec.getStart()/60f, rangeInSec.getStop()/60f}, new double[] {0.0, maxY, maxY}, GraphType.area, "Boundaries", new Color(102, 204, 255, 50), 4.0f));
		}
		
		return traces;
	}

	private ArrayList<XYTrace> extractFragmentTraces(FragmentIon[] primaryIonObjects, String sampleName, QuantitativeDIAData quantitativeData, StripeFileInterface file) throws IOException, SQLException {
		Pair<ArrayList<XYTrace>, Double> tracesAndMaxY;
		ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
		TObjectDoubleHashMap<FragmentIon> targetIonObjects=new TObjectDoubleHashMap<FragmentIon>();
		ArrayList<FragmentIon> offtargetIonObjects=new ArrayList<FragmentIon>();
		
		XYTrace quantitativePeaks=new XYTrace(quantitativeData.getMassArray(), quantitativeData.getIntensityArray(), GraphType.spectrum, sampleName);
		Collections.sort(quantitativePeaks.getPoints());
		Pair<double[], double[]> peaksArrays=quantitativePeaks.toArrays();
		double[] targets=peaksArrays.x;
		double[] intensities=peaksArrays.y;
		
		for (FragmentIon ion : primaryIonObjects) {
			Optional<Integer> index=parameters.getFragmentTolerance().getIndex(targets, ion.getMass());
			if (index.isPresent()) {
				targetIonObjects.put(ion, intensities[index.get()]);
			} else {
				offtargetIonObjects.add(ion);
			}
		}
		FragmentIon[] targetIonArray=targetIonObjects.keys(new FragmentIon[targetIonObjects.size()]);
		FragmentIon[] offTargetIonArray=offtargetIonObjects.toArray(new FragmentIon[offtargetIonObjects.size()]);
		
		Range rangeInSec=quantitativeData.getRtScanRange();
		Range rangeInMins=new Range(rangeInSec.getStart()/60f, rangeInSec.getStop()/60f);
		ArrayList<FragmentScan> stripes=file.getStripes(quantitativeData.getPrecursorMZ(), rangeInSec.getStart()-RT_EXTRACTION_MARGIN_IN_SEC, rangeInSec.getStop()+RT_EXTRACTION_MARGIN_IN_SEC, false);
		
		ArrayList<Spectrum> downcastedSpectra=FragmentScan.downcastStripeToSpectrum(stripes);

		HashMap<FragmentIon, XYTrace> targetFragmentTraceMap=ChromatogramExtractor.extractFragmentChromatograms(parameters.getFragmentTolerance(), targetIonArray, downcastedSpectra, null,
				GraphType.boldline, sgSmoothBox.isSelected(), backgroundSubtractBox.isSelected());
		HashMap<FragmentIon, XYTrace> offTargetFragmentTraceMap=ChromatogramExtractor.extractFragmentChromatograms(parameters.getFragmentTolerance(), offTargetIonArray, downcastedSpectra,
				null, GraphType.dashedline, sgSmoothBox.isSelected(), backgroundSubtractBox.isSelected());

		traces.addAll(targetFragmentTraceMap.values());
		double maxY=0.1;
		
		for (Entry<FragmentIon, XYTrace> ionEntry : targetFragmentTraceMap.entrySet()) {
			XYTrace trace=ionEntry.getValue();
			XYPoint xy=trace.getMaxXYInRange(rangeInMins);
			if (xy!=null) {
				if (xy.getY()>maxY) {
					maxY=xy.getY();
				}
				double intensity=targetIonObjects.get(ionEntry.getKey());
				//traces.add(new XYTrace(new double[] {xy.x}, new double[] {xy.y}, GraphType.text, trace.getName()+" ("+formatter.format(intensity).toLowerCase()+")"));
			}
		}
		
		traces.addAll(offTargetFragmentTraceMap.values());
		
		if (traces.size()>0) {
			// extra 0 point in case there is no data shown (or all 0s)
			traces.add(new XYTrace(new double[] {rangeInMins.getStart()-Float.MIN_VALUE, rangeInMins.getStart(), rangeInMins.getStop()}, new double[] {0.0, maxY, maxY},
					GraphType.area, "Boundaries", new Color(102, 204, 255, 50), 4.0f));
		}
		
		return traces;
	}

}
