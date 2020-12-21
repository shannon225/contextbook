package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.text.AttributedString;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.zip.DataFormatException;
import java.util.Optional;

import javax.swing.BoxLayout;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.RowFilter;
import javax.swing.SwingUtilities;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.CombinedRangeXYPlot;
import org.jfree.ui.TextAnchor;
import org.relaxng.datatype.DatatypeException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.LibraryReportExtractor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideReportData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LabeledComponent;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingWorkerProgress;
import edu.washington.gs.maccoss.encyclopedia.gui.massspec.ChromatogramCharter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.StringUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.AcquiredSpectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.ChromatogramExtractor;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.QuantitativeDIAData;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.map.hash.TObjectDoubleHashMap;

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
	private final JCheckBox plotPrecursors=new JCheckBox("Plot Precursors");
	
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
		plotPrecursors.setBackground(LabeledComponent.BACKGROUND_COLOR);
		checkBoxes.add(plotPrecursors);
		options.add(checkBoxes);
		simplifyPlots.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				updateToSelectedPeptide();
			}
		});
		plotPrecursors.addActionListener(new ActionListener() {
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
		if (library!=null) {
			try {
			ArrayList<LibraryEntry> entries=library.getEntries(peptide.getPeptideModSeq(), peptide.getPrecursorCharge(), false);
			ArrayList<Spectrum> newDataArray=new ArrayList<>();
			newDataArray.addAll(entries);
			newDataArray.addAll(Arrays.asList(dataArray));
			dataArray=newDataArray.toArray(new Spectrum[newDataArray.size()]);
			
			ArrayList<String> newCategories=new ArrayList<>();
			if (entries.size()==1) {
				newCategories.add("Library");
			} else if (entries.size()>1) {
				for (int i = 1; i <= entries.size(); i++) {
					newCategories.add("Library "+i);
				}
			}
			newCategories.addAll(Arrays.asList(categories));
			categories=newCategories.toArray(new String[newCategories.size()]);
			} catch (Exception e) {
				Logger.errorLine("Error reading library entry!");
				Logger.errorException(e);
			}
		}
		
		return FragmentIonConsistencyCharter.getBarChart(peptide, dataArray, categories, parameters);
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

				System.out.println(split1.getDividerLocation()+", "+split2.getDividerLocation());
				if (split1.getDividerLocation()<=40||split2.getDividerLocation()<=40) {
					split2.setDividerLocation(200);
					split1.setDividerLocation(200);
				}
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
			System.out.println(origSampleNames[i]+" --> "+dataArray[i]);
		}
		String[] sampleNames=StringUtils.getUniquePortion(origSampleNames);
		barChart.setChart(getBarChart(sampleNames, totalTICs).getChart());
		
		stackedBarChart.setChart(getStackedBarChart(entry, sampleNames, dataArray).getChart());
		
		int cols=(Integer)numberOfColumns.getSelectedItem();
		boolean simplify=simplifyPlots.isSelected();
		if (files.size()<cols) cols=files.size();
		
		JPanel right=new JPanel(new GridLayout(0, simplify?1:cols));
		right.setBackground(Color.WHITE);
		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		FragmentIon[] primaryIonObjects=model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), true);
		Logger.logLine("Graphing "+entry.getPeptideModSeq()+" ("+primaryIonObjects.length+")"+"...");
		
		try {
			
			double globalMaxYFragment=0.0;
			double globalMaxYPrecursor=0.0;
			
			ArrayList<ArrayList<XYTrace>> allFragmentTraces=new ArrayList<ArrayList<XYTrace>>();
			ArrayList<ArrayList<XYTrace>> allPrecursorTraces=new ArrayList<ArrayList<XYTrace>>();
			for (int i=0; i<origSampleNames.length; i++) {
				String sampleName = origSampleNames[i];
				QuantitativeDIAData quantitativeData=entry.getQuantitativeData(sampleName);
				StripeFileInterface file=files.get(i);

				ArrayList<XYTrace> fragmentTraces;
				ArrayList<XYTrace> precursorTraces;
				if (quantitativeData==null) {
					fragmentTraces=null;
					precursorTraces=null;
				} else {
					precursorTraces=extractPrecursorTraces(sampleName, quantitativeData, file);
					fragmentTraces=extractFragmentTraces(primaryIonObjects, sampleName, quantitativeData, file);
				}

				for (XYTrace xyTrace : fragmentTraces) {
					if (xyTrace.getType()==GraphType.line) {
						globalMaxYFragment=Math.max(globalMaxYFragment, xyTrace.getMaxY());
					}
				}

				for (XYTrace xyTrace : precursorTraces) {
					if (xyTrace.getType()==GraphType.line) {
						globalMaxYPrecursor=Math.max(globalMaxYPrecursor, xyTrace.getMaxY());
					}
				}
				allFragmentTraces.add(fragmentTraces);
				allPrecursorTraces.add(precursorTraces);
			}

			globalMaxYFragment=globalMaxYFragment*1.05;
			globalMaxYPrecursor=globalMaxYPrecursor*1.05;

			CombinedRangeXYPlot parent=null;
			for (int i=0; i<sampleNames.length; i++) {
				if (simplify) {
					ChartPanel fragmentChart=Charter.getChart("Retention Time (min)", "Intensity", false, globalMaxYFragment, allFragmentTraces.get(i).toArray(new XYTrace[0]));
					
					fragmentChart.getChart().getXYPlot().clearAnnotations();
					
					ValueAxis domainAxis = fragmentChart.getChart().getXYPlot().getDomainAxis();
					ValueAxis rangeAxis = fragmentChart.getChart().getXYPlot().getRangeAxis();
					XYTextAnnotation annotation = new XYTextAnnotation(sampleNames[i], domainAxis.getLowerBound(), rangeAxis.getUpperBound());
					annotation.setTextAnchor(TextAnchor.TOP_LEFT);
					annotation.setFont(new Font("News Gothic MT", Font.BOLD, 14));
					fragmentChart.getChart().getXYPlot().addAnnotation(annotation);
					
					domainAxis.setLabel(null);
					domainAxis.setTickLabelFont(new Font("News Gothic MT", Font.PLAIN, 12));
					
					if (i%cols==0) {
						// ADD label the domain of the left most plots (FIXME FIND A BETTER WAY TO DO THIS)
				        parent = new CombinedRangeXYPlot(rangeAxis);

						final ChartPanel chartPanel=new ChartPanel(new JFreeChart(parent), false);
						chartPanel.getChart().removeLegend();
						chartPanel.getChart().setBackgroundPaint(Color.white);
						chartPanel.setMinimumDrawWidth(0);
						chartPanel.setMinimumDrawHeight(0);
						chartPanel.setMaximumDrawWidth(Integer.MAX_VALUE);
						chartPanel.setMaximumDrawHeight(Integer.MAX_VALUE);
				        right.add(chartPanel);
				        
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
					ChartPanel panel=ChromatogramCharter.createChart(Optional.ofNullable(allPrecursorTraces.get(i)),
							Optional.ofNullable(allFragmentTraces.get(i)), globalMaxYPrecursor, 
							globalMaxYFragment);
					panel.getChart().setTitle(sampleNames[i]);
					right.add(panel);
				}
			}
	
			split.setRightComponent(right);

		} catch (DataFormatException sqle) {
			Logger.errorLine("Error reading raw files!");
			Logger.errorException(sqle);
		} catch (DatatypeException sqle) {
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
	
	private ArrayList<XYTrace> extractPrecursorTraces(String sampleName, QuantitativeDIAData quantitativeData, StripeFileInterface file) throws IOException, SQLException, DatatypeException, DataFormatException {
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
		XYTraceInterface[] traceArray=ChromatogramExtractor.extractPrecursorChromatograms(parameters.getPrecursorTolerance(), quantitativeData.getPrecursorMZ(), quantitativeData.getPrecursorCharge(), precursors);

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
				GraphType.boldline);
		HashMap<FragmentIon, XYTrace> offTargetFragmentTraceMap=ChromatogramExtractor.extractFragmentChromatograms(parameters.getFragmentTolerance(), offTargetIonArray, downcastedSpectra,
				null, GraphType.dashedline);

		traces.addAll(targetFragmentTraceMap.values());
		double maxY=0.0;
		
		for (Entry<FragmentIon, XYTrace> ionEntry : targetFragmentTraceMap.entrySet()) {
			XYTrace trace=ionEntry.getValue();
			XYPoint xy=trace.getMaxXYInRange(rangeInMins);
			if (xy.getY()>maxY) {
				maxY=xy.getY();
			}
			double intensity=targetIonObjects.get(ionEntry.getKey());
			traces.add(new XYTrace(new double[] {xy.x}, new double[] {xy.y}, GraphType.text, trace.getName()+" ("+formatter.format(intensity).toLowerCase()+")"));
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
