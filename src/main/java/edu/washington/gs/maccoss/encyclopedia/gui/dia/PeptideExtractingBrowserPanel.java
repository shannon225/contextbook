package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.JToolBar;
import javax.swing.RowFilter;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.UIManager;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.event.PlotChangeEvent;
import org.jfree.chart.event.PlotChangeListener;
import org.jfree.chart.plot.XYPlot;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ExpectedFragmentationScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneFragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanRawScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScoringTask;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SimplePeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.ExtendedChartPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.massspec.ChromatogramCharter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.AcquiredSpectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.ChromatogramExtractor;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SpectrumUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class PeptideExtractingBrowserPanel extends JPanel {
	private static final long serialVersionUID=1L;
	
	private final PecanRawScorer scorer;

	private final SearchParameters parameters;
	private final FileChooserPanel diaFile;
	private final JTextField peptide=new JTextField("VATVSLPR");
	private final SpinnerModel charge=new SpinnerNumberModel(2, 1, 5, 1);
	private final JSplitPane chromatogramSplit=new JSplitPane(JSplitPane.VERTICAL_SPLIT);
	private final JSplitPane spectrumSplit=new JSplitPane(JSplitPane.VERTICAL_SPLIT);
	private final JSplitPane horizontalSplit=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	private final JTable table;
	private final TableRowSorter<TableModel> rowSorter;
	private final JTextField jtfFilter;
	private final DIAScanTableModel model;
	private final JCheckBox sgSmoothBox;
	private final JCheckBox backgroundSubtractBox;

	private StripeFileInterface dia=null;

	public static void main(String[] args) {
		try {
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		} catch (Exception e) {
			Logger.errorLine("Error setting look and feel!");
			Logger.errorException(e);
		}
		OS os=OSDetector.getOS();
		switch (os) {
			case MAC:
				System.setProperty("com.apple.mrj.application.apple.menu.about.name", "DIA Browser");
				System.setProperty("apple.laf.useScreenMenuBar", "true");
				break;

			default:
				break;
		}

		EventQueue.invokeLater(new Runnable() {
			@Override
			public void run() {
				final JFrame f=new JFrame("Peptide Browser");
				f.addWindowListener(new WindowAdapter() {
					public void windowClosing(WindowEvent e) {
						System.exit(0);
					}
				});

				SearchParameters params=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), DataAcquisitionType.OVERLAPPING_DIA, false, true, false);
				f.getContentPane().add(new PeptideExtractingBrowserPanel(params), BorderLayout.CENTER);

				f.pack();
				f.setSize(new Dimension(1900, 1030)); // for 1920x1080
				f.setVisible(true);
			}
		});

		Logger.logLine("Launching DIA Browser");
	}
	
	public PeptideExtractingBrowserPanel(SearchParameters parameters) {
		super(new BorderLayout());
		this.parameters=parameters; 
		scorer=new PecanRawScorer(parameters.getFragmentTolerance(), new ExpectedFragmentationScorer(parameters, 3));

		diaFile=new FileChooserPanel(null, "RAW File", StripeFileGenerator.getFilenameFilter(), true) {
			private static final long serialVersionUID=1L;

			@Override
			public void update(File... filename) {
				super.update(filename);
				if (filename!=null&&filename.length>0&&filename[0]!=null) {
					try {
						Logger.logLine("Reading file...");

						dia=StripeFileGenerator.getFile(filename[0], PeptideExtractingBrowserPanel.this.parameters);
						
						Logger.logLine("Finished reading file.");
						resetPeptide(peptide.getText(), (Integer) charge.getValue());
					} catch (Exception e) {
						e.printStackTrace();
						JOptionPane.showMessageDialog(PeptideExtractingBrowserPanel.this, "Sorry, there was a problem reading ["+filename[0].getName()+"]: "+e.getMessage(), "Error Opening DIA File",
								JOptionPane.ERROR_MESSAGE);
					}
				} else {
					dia=null;
				}
			}
		};

		sgSmoothBox=new JCheckBox("Smooth");
		sgSmoothBox.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				resetPeptide(peptide.getText(), (Integer) charge.getValue());
			}
		});
		backgroundSubtractBox=new JCheckBox("Background Subtract");
		backgroundSubtractBox.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				resetPeptide(peptide.getText(), (Integer) charge.getValue());
			}
		});

		JToolBar bar=new JToolBar();
		bar.add(diaFile);
		bar.add(sgSmoothBox);
		bar.add(backgroundSubtractBox);
		bar.add(new JLabel("Peptide Sequence:"));
		bar.add(peptide);
		bar.add(new JLabel("Charge:"));
		bar.add(new JSpinner(charge));

		JButton button=new JButton("GO");
		button.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				resetPeptide(peptide.getText(), (Integer) charge.getValue());
			}
		});
		bar.add(button);
		bar.add(new JPanel());

		add(bar, BorderLayout.NORTH);
		
		model=new DIAScanTableModel();
		table=new JTable(model) {
			private static final long serialVersionUID=1L;

			@Override
			public Object getValueAt(int row, int column) {
				if (column==0) return row+1;
				return super.getValueAt(row, column);
			}
		};
		rowSorter=new TableRowSorter<TableModel>(table.getModel());
		table.setRowSorter(rowSorter);
		
		table.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			@Override
			public void valueChanged(ListSelectionEvent e) {
				updateToSelected();
			}
		});

		jtfFilter=new JTextField();
		jtfFilter.getDocument().addDocumentListener(new DocumentListener() {
			@Override
			public void insertUpdate(DocumentEvent e) {
				String text=jtfFilter.getText();

				System.out.println("FILTER: "+text);
				if (text.trim().length()==0) {
					rowSorter.setRowFilter(null);
				} else {
					rowSorter.setRowFilter(RowFilter.regexFilter("(?i)"+text));
				}
			}

			@Override
			public void removeUpdate(DocumentEvent e) {
				String text=jtfFilter.getText();

				if (text.trim().length()==0) {
					rowSorter.setRowFilter(null);
				} else {
					rowSorter.setRowFilter(RowFilter.regexFilter("(?i)"+text));
				}
			}

			@Override
			public void changedUpdate(DocumentEvent e) {
				throw new UnsupportedOperationException("Not supported yet.");
			}
		});


		JPanel searchPanel=new JPanel(new BorderLayout());
		searchPanel.add(new JLabel("Search:"), BorderLayout.WEST);
		searchPanel.add(jtfFilter, BorderLayout.CENTER);
		
		spectrumSplit.setTopComponent(new JScrollPane(table));
		
		horizontalSplit.setLeftComponent(chromatogramSplit);
		horizontalSplit.setRightComponent(spectrumSplit);
		add(horizontalSplit, BorderLayout.CENTER);
	}


	public void updateToSelected() {
		int[] selection=table.getSelectedRows();
		if (selection.length<=0) return;
		
		ArrayList<AcquiredSpectrum> entries=new ArrayList<AcquiredSpectrum>();
		for (int row : selection) {
			AcquiredSpectrum entry=model.getSelectedRow(table.convertRowIndexToModel(row));
			entries.add(entry);
		}
		resetScan(entries);
	}

	public void resetScan(ArrayList<AcquiredSpectrum> entries) {
		int locationSpectrum=spectrumSplit.getDividerLocation();
		if (locationSpectrum<=10) {
			locationSpectrum=getHeight()/2;
		} else if (locationSpectrum>=getHeight()*.9f) {
			locationSpectrum=getHeight()/2;
		}
		
		final String peptideModSeq=peptide.getText();
		final byte precursorCharge=((Integer)charge.getValue()).byteValue();
		
		final Spectrum spectrum;
		if (entries.size()==1) {
			spectrum=entries.get(0);
		} else {
			spectrum=SpectrumUtils.mergeSpectra(downcast(entries), parameters.getFragmentTolerance());
		}
		
		SimplePeptidePrecursor precursor=new SimplePeptidePrecursor(peptideModSeq, precursorCharge, parameters.getAAConstants());
		AnnotatedLibraryEntry entry=new AnnotatedLibraryEntry(precursor, spectrum, parameters);

		final ChartPanel spectrumChart=Charter.getChart(entry);
		spectrumSplit.setBottomComponent(spectrumChart);

		spectrumSplit.setDividerLocation(locationSpectrum);
	}
	
	private ArrayList<Spectrum> downcast(ArrayList<AcquiredSpectrum> spectra) {
		ArrayList<Spectrum> ret=new ArrayList<>();
		ret.addAll(spectra);
		return ret;
	}
	
	volatile double lowerBound=0.0;
	volatile double upperBound=0.0;

	public void resetPeptide(String peptide, int charge) {
		int location=horizontalSplit.getDividerLocation();
		if (location<=10) {
			location=getWidth()/2;
		}
		int locationChromatogram=chromatogramSplit.getDividerLocation();
		if (locationChromatogram<=10) {
			locationChromatogram=getHeight()/2;
		}
		int locationSpectrum=spectrumSplit.getDividerLocation();
		if (locationSpectrum<=10) {
			locationSpectrum=getHeight()/2;
		}
		
		if (peptide==null||peptide.length()==0||dia==null) {
			chromatogramSplit.setTopComponent(new JLabel("Select a peptide!"));
			chromatogramSplit.setBottomComponent(new JPanel());
			spectrumSplit.setTopComponent(new JPanel());
			spectrumSplit.setBottomComponent(new JPanel());
		} else {
			Logger.logLine("Parsing peptide...");
			PecanOneFragmentationModel fragModel=new PecanOneFragmentationModel(new FastaPeptideEntry(peptide), parameters.getAAConstants());
			ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
			AnnotatedLibraryEntry entry=fragModel.getUnitSpectrum(dia.getOriginalFileName(), new HashSet<String>(), (byte)charge, 0.0f, parameters);
			entries.add(entry);
			
			try {
				ArrayList<PrecursorScan> precursors=dia.getPrecursors(0.0f, Float.MAX_VALUE);
				ArrayList<FragmentScan> stripes=dia.getStripes(entry.getPrecursorMZ(), 0.0f, Float.MAX_VALUE, false);
				
				ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
				ArrayList<XYTrace> precursorTraces=new ArrayList<XYTrace>();
//				for (Entry<LibraryEntry, PeptideScoringResult> resultEntry : result.entrySet()) {
//					FragmentationScoringResult peptideResult=(FragmentationScoringResult)resultEntry.getValue();
//
//					for (XYTrace trace : peptideResult.getFragmentationTraces()) {
//						//XYTrace sgSmoothed=SkylineSGFilter.paddedSavitzkyGolaySmooth(trace);
//						traces.add(trace);
//					}
//				}

				ArrayList<AcquiredSpectrum> scans=new ArrayList<AcquiredSpectrum>();
				for (FragmentScan scan : stripes) {
					scans.add(scan);
				}
				model.updateEntries(scans);
				if (scans.size()>0) {
					table.addRowSelectionInterval(0, 0);
				} else {
				}
				updateToSelected();

				HashMap<FragmentIon, XYTrace> targetFragmentTraceMap=ChromatogramExtractor.extractFragmentChromatograms(parameters.getFragmentTolerance(), entry.getIonAnnotations(), stripes, null,
						GraphType.boldline, sgSmoothBox.isSelected(), backgroundSubtractBox.isSelected());
				traces.addAll(targetFragmentTraceMap.values());

				XYTraceInterface[] precursorTraceArray=ChromatogramExtractor.extractPrecursorChromatograms(parameters.getPrecursorTolerance(), 
						entry.getPrecursorMZ(), entry.getPrecursorCharge(), precursors, sgSmoothBox.isSelected(), backgroundSubtractBox.isSelected());
				
				double globalMaxYPrecursor=0.0;
				for (XYTraceInterface trace : precursorTraceArray) {
					if (trace instanceof XYTrace) {
						globalMaxYPrecursor=Math.max(globalMaxYPrecursor, General.max(trace.toArrays().y));
						precursorTraces.add((XYTrace)trace);
					}
				}

				double globalMaxYFragment=0.0;
				for (XYTrace trace : traces) {
					globalMaxYFragment=Math.max(globalMaxYFragment, General.max(trace.toArrays().y));
				}

				ExtendedChartPanel chart=ChromatogramCharter.createChart(Optional.ofNullable(precursorTraces),
						Optional.ofNullable(traces), globalMaxYPrecursor, globalMaxYFragment);
				//ExtendedChartPanel chart=Charter.getChart("Retention Time (min)", "Intensity", false, traces.toArray(new XYTrace[traces.size()]));
				addAnnotations(targetFragmentTraceMap, chart);
				chart.getChart().getXYPlot().addChangeListener(new PlotChangeListener() {
					
					@Override
					public void plotChanged(PlotChangeEvent event) {
						
						synchronized (this) {
							addAnnotations(targetFragmentTraceMap, chart);
						}
					}
				});
				chromatogramSplit.setTopComponent(chart);
				

				
				/*BlockingQueue<PeptideScoringResult> ionCountResultsQueue=new LinkedBlockingQueue<PeptideScoringResult>();
				IonCountingScoringTask ionCount=new IonCountingScoringTask(scorer, entries, stripes, 2.5f, new PrecursorScanMap(new ArrayList<PrecursorScan>()), ionCountResultsQueue, parameters);
				ionCount.call();
				
				PeptideScoringResult ionCountResult=ionCountResultsQueue.take();
				XYTraceInterface xytrace=ionCountResult.getTrace();
				Pair<double[], double[]> trace=xytrace.toArrays();
				double[] newx=General.multiply(trace.x, 1.0f/60.0f); // scale to minutes
				XYTraceInterface ionCounttrace=new XYTrace(newx, trace.y, xytrace.getType(), xytrace.getName(), xytrace.getColor(), xytrace.getThickness());
				
				ChartPanel ionCountchart=Charter.getChart("Retention Time (min)", "RawScore", false, ionCounttrace);*/
				BlockingQueue<AbstractScoringResult> resultsQueue=new LinkedBlockingQueue<AbstractScoringResult>();
				XCorDIAOneScoringTask xcorrTask=new XCorDIAOneScoringTask(new XCorDIAOneScorer(parameters, null), entries, stripes, new Range(0.0f, 0.0f), 2.5f, new PrecursorScanMap(new ArrayList<PrecursorScan>()), resultsQueue, parameters);
				xcorrTask.call();
				
				AbstractScoringResult ionCountResult=resultsQueue.take();
				XYTraceInterface xytrace=ionCountResult.getTrace();
				Pair<double[], double[]> trace=xytrace.toArrays();
				double[] newx=General.multiply(trace.x, 1.0f/60.0f); // scale to minutes
				XYTraceInterface ionCounttrace=new XYTrace(newx, trace.y, xytrace.getType(), xytrace.getName(), xytrace.getColor(), xytrace.getThickness());
				
				ChartPanel ionCountchart=Charter.getChart("Retention Time (min)", "XCorr", false, ionCounttrace);

				chromatogramSplit.setBottomComponent(ionCountchart);
				
				
			} catch (Exception e) {
				e.printStackTrace();
				JOptionPane.showMessageDialog(PeptideExtractingBrowserPanel.this, "Sorry, there was a problem reading the precursor window that contains ["+entry.getPrecursorMZ()+"]: "+e.getMessage(), "Error Reading DIA File",
						JOptionPane.ERROR_MESSAGE);
			}
			Logger.logLine("Finished reading peptide "+peptide);
		}

		spectrumSplit.setDividerLocation(locationSpectrum);
		chromatogramSplit.setDividerLocation(locationChromatogram);
		horizontalSplit.setDividerLocation(location);
	}

	private void addAnnotations(HashMap<FragmentIon, XYTrace> targetFragmentTraceMap, ExtendedChartPanel chart) {
		XYPlot plot = chart.getChart().getXYPlot();
		org.jfree.data.Range jfreeRange=plot.getDomainAxis().getRange();
		if (jfreeRange.getLowerBound()!=lowerBound||jfreeRange.getUpperBound()!=upperBound) {
			lowerBound=jfreeRange.getLowerBound();
			upperBound=jfreeRange.getUpperBound();
			System.out.println("Found: "+lowerBound+" / "+upperBound);
			
			Range range=new Range(jfreeRange.getLowerBound(), jfreeRange.getUpperBound());
			for (Entry<FragmentIon, XYTrace> ionEntry : targetFragmentTraceMap.entrySet()) {
				FragmentIon ion=ionEntry.getKey();
				XYTrace trace=ionEntry.getValue();
				XYPoint xy=trace.getMaxXYInRange(range);

				if (xy!=null) {
					XYTextAnnotation annotation=new XYTextAnnotation(ion.getName()+" ("+String.format("%.1f", ion.getMass())+" m/z)", xy.x, xy.y/chart.getDivider()*1.01);
					System.out.println("      "+trace.getName()+" = "+xy.x+" / "+(xy.y/chart.getDivider()));
					plot.addAnnotation(annotation);
				}
			}
		}
	}
}
