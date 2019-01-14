package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Frame;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Optional;

import javax.swing.BoxLayout;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
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

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.LibraryReportExtractor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideReportData;
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
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.StringUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.ChromatogramExtractor;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.QuantitativeDIAData;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import gnu.trove.map.hash.TObjectDoubleHashMap;

public class MultiResultsBrowserPanel extends JPanel {
	private static final long serialVersionUID=1L;

	private static final int RT_EXTRACTION_MARGIN_IN_SEC=15;
	
	private final FileChooserPanel elibFileChooser;
	private final JSplitPane split=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	
	private final JTable sampleTable;
	private final TableRowSorter<TableModel> sampleRowSorter;
	private final SampleTableModel sampleModel;
	private final JTable peptideTable;
	private final TableRowSorter<TableModel> peptideRowSorter;
	private final JTextField jtfFilter;
	private final MultiPeptideResultsTableModel peptideModel;
	private final SearchParameters parameters;
	private final ChartPanel barChart;
	private final JComboBox<Integer> minimumNumberOfFragments=new JComboBox<Integer>(new Integer[] {0, 1, 2, 3, 4, 5});
	private final JComboBox<Integer> numberOfColumns=new JComboBox<Integer>(new Integer[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
	
	private final int defaultMinimumNumberOfFragmentsIndex=3;
	private final int defaultNumberOfColumnsIndex=1;
	
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
			private static final long serialVersionUID=1L;

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
				String text=jtfFilter.getText();

				if (text.trim().length()==0) {
					peptideRowSorter.setRowFilter(null);
				} else {
					peptideRowSorter.setRowFilter(RowFilter.regexFilter("(?i)"+text));
				}
			}

			@Override
			public void removeUpdate(DocumentEvent e) {
				String text=jtfFilter.getText();

				if (text.trim().length()==0) {
					peptideRowSorter.setRowFilter(null);
				} else {
					peptideRowSorter.setRowFilter(RowFilter.regexFilter("(?i)"+text));
				}
			}

			@Override
			public void changedUpdate(DocumentEvent e) {
				throw new UnsupportedOperationException("Not supported yet.");
			}
		});

		peptideTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			@Override
			public void valueChanged(ListSelectionEvent e) {
				updateToSelectedPeptide();
			}
		});

		elibFileChooser=new FileChooserPanel(null, "Library", new SimpleFilenameFilter(LibraryFile.DLIB, LibraryFile.ELIB), true) {
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
		options.add(new LabeledComponent("Minimum # Fragments", minimumNumberOfFragments));
		options.add(new LabeledComponent("Number of Columns", numberOfColumns));
		
		JPanel tablePanel=new JPanel(new GridLayout(0, 1));
		tablePanel.add(new JScrollPane(sampleTable, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED));
		tablePanel.add(new JScrollPane(peptideTable, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED));
		barChart=getBarChart(new String[] {}, new float[] {});
		tablePanel.add(barChart, BorderLayout.SOUTH);

		JPanel searchPanel=new JPanel(new BorderLayout());
		searchPanel.add(new JLabel("Search:"), BorderLayout.WEST);
		searchPanel.add(jtfFilter, BorderLayout.CENTER);

		JPanel left=new JPanel(new BorderLayout());
		left.add(options, BorderLayout.NORTH);
		left.add(tablePanel, BorderLayout.CENTER);
		left.add(searchPanel, BorderLayout.SOUTH);
		
		
		
		split.setLeftComponent(left);
		split.setRightComponent(new JLabel("Select a peptide!"));
		

		setLayout(new BorderLayout());
		add(split, BorderLayout.CENTER);
	}

	public ChartPanel getBarChart(String[] categories, float[] intensities) {
		return Charter.getBarChart(null, "Sample", "Total Intensity", categories, intensities);
	}

	public void askForLibrary() {
		elibFileChooser.askForFiles();
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
			}
		};
		worker.execute();
	}
	
	public void updateTables(final File f, final Pair<ArrayList<String>, ArrayList<PeptideReportData>> pair) {
		SwingWorkerProgress<Pair<ArrayList<StripeFileInterface>, ArrayList<PeptideReportData>>> worker=new SwingWorkerProgress<Pair<ArrayList<StripeFileInterface>, ArrayList<PeptideReportData>>>((Frame)SwingUtilities.getWindowAncestor(this), "Please wait...", "Reading Library") {
			@Override
			protected Pair<ArrayList<StripeFileInterface>, ArrayList<PeptideReportData>> doInBackgroundForReal() throws Exception {
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
		for (int i=0; i<origSampleNames.length; i++) {
			origSampleNames[i]=files.get(i).getOriginalFileName();
			QuantitativeDIAData quantitativeData=entry.getQuantitativeData(origSampleNames[i]);
			if (quantitativeData!=null) {
				totalTICs[i]=quantitativeData.getTIC();
			}
		}
		String[] sampleNames=StringUtils.getUniquePortion(origSampleNames);
		barChart.setChart(getBarChart(sampleNames, totalTICs).getChart());
		int cols=(Integer)numberOfColumns.getSelectedItem();
		if (files.size()<cols) cols=files.size();
		
		JPanel right=new JPanel(new GridLayout(0, cols));
		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		double precursorMz=parameters.getAAConstants().getChargedMass(entry.getPeptideModSeq(), entry.getPrecursorCharge());
		
		try {
			FragmentIon[] primaryIonObjects=model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), true);
			Logger.logLine("Graphing "+entry.getPeptideModSeq()+" ("+primaryIonObjects.length+")"+"...");
			
			double globalMaxY=0.0;
			ArrayList<ChartPanel> allPanels=new ArrayList<ChartPanel>();
			
			ArrayList<ArrayList<XYTrace>> allTraces=new ArrayList<ArrayList<XYTrace>>();
			for (int i=0; i<origSampleNames.length; i++) {
				ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
				QuantitativeDIAData quantitativeData=entry.getQuantitativeData(origSampleNames[i]);
				
				if (quantitativeData!=null) {
					TObjectDoubleHashMap<FragmentIon> targetIonObjects=new TObjectDoubleHashMap<FragmentIon>();
					ArrayList<FragmentIon> offtargetIonObjects=new ArrayList<FragmentIon>();
					
					for (int j=0; j<quantitativeData.getMassArray().length; j++) {
						System.out.println(quantitativeData.getMassArray()[j]+"\t"+quantitativeData.getIntensityArray()[j]);
					}
					System.out.println();
					XYTrace quantitativePeaks=new XYTrace(quantitativeData.getMassArray(), quantitativeData.getIntensityArray(), GraphType.spectrum, origSampleNames[i]);
					Collections.sort(quantitativePeaks.getPoints());
					Pair<double[], double[]> peaksArrays=quantitativePeaks.toArrays();
					double[] targets=peaksArrays.x;
					double[] intensities=peaksArrays.y;
					
					for (FragmentIon ion : primaryIonObjects) {
						Optional<Integer> index=parameters.getFragmentTolerance().getIndex(targets, ion.mass);
						if (index.isPresent()) {
							targetIonObjects.put(ion, intensities[index.get()]);
						} else {
							offtargetIonObjects.add(ion);
						}
					}
					FragmentIon[] targetIonArray=targetIonObjects.keys(new FragmentIon[targetIonObjects.size()]);
					FragmentIon[] offTargetIonArray=offtargetIonObjects.toArray(new FragmentIon[offtargetIonObjects.size()]);
					
					StripeFileInterface file=files.get(i);
					Range rangeInSec=quantitativeData.getRtScanRange();
					Range rangeInMins=new Range(rangeInSec.getStart()/60f, rangeInSec.getStop()/60f);
					ArrayList<FragmentScan> stripes=file.getStripes(precursorMz, rangeInSec.getStart()-RT_EXTRACTION_MARGIN_IN_SEC, rangeInSec.getStop()+RT_EXTRACTION_MARGIN_IN_SEC, false);
					
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
					
					globalMaxY=Math.max(globalMaxY, maxY);
					traces.addAll(offTargetFragmentTraceMap.values());
					
					if (traces.size()>0) {
						// extra 0 point in case there is no data shown (or all 0s)
						traces.add(new XYTrace(new double[] {rangeInMins.getStart()-Float.MIN_VALUE, rangeInMins.getStart(), rangeInMins.getStop()}, new double[] {0.0, maxY, maxY},
								GraphType.area, "Boundaries", new Color(102, 204, 255, 50), 4.0f));
					}
				}
				allTraces.add(traces);
			}
	
			if (globalMaxY>0.0) {
				globalMaxY=globalMaxY*1.05;
			}

			for (int i=0; i<sampleNames.length; i++) {
				ArrayList<XYTrace> traces=allTraces.get(i);
				ChartPanel fragmentChart=Charter.getChart("Retention Time (min)", "Intensity", false, globalMaxY, traces.toArray(new XYTrace[traces.size()]));
				fragmentChart.getChart().setTitle(sampleNames[i]);
				allPanels.add(fragmentChart);
				right.add(fragmentChart);
			}
	
			split.setRightComponent(right);
			
		} catch (SQLException sqle) {
			Logger.errorLine("Error reading raw files!");
			Logger.errorException(sqle);
		} catch (IOException ioe) {
			Logger.errorLine("Error reading raw files!");
			Logger.errorException(ioe);
		}
		
		split.setDividerLocation(location);
	}
}
