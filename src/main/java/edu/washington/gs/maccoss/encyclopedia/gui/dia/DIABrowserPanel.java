package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.Frame;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Optional;

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.RowFilter;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LabeledComponent;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingWorkerProgress;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.AcquiredSpectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.ChromatogramExtractor;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Polymer;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PolymerIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SpectrumComparator;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SpectrumUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.PivotTableGenerator;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TDoubleDoubleHashMap;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class DIABrowserPanel extends JPanel {
	private static final String STRUCTURE_TITLE="Structure";
	private static final String INTENSITY_DISTRIBUTION_TITLE="Intensity Distributions";
	private static final String BOXPLOT_TITLE="Range Statistics";
	private static final long serialVersionUID=1L;
	public static final Color[] colors=new Color[] {Color.red, Color.blue, Color.green, Color.cyan, Color.magenta, Color.orange, Color.yellow, Color.pink, Color.gray, 
			Color.red.darker(), Color.blue.darker(), Color.green.darker(), Color.cyan.darker(), Color.magenta.darker(), Color.orange.darker(), Color.yellow.darker(), Color.pink.darker(), Color.gray.darker()};

	private final FileChooserPanel rawFileChooser;
	private final JSplitPane boxplotSplit=new JSplitPane(JSplitPane.VERTICAL_SPLIT);
	private final JSplitPane distributionSplit=new JSplitPane(JSplitPane.VERTICAL_SPLIT);
	private final JSplitPane rawSplit=new JSplitPane(JSplitPane.VERTICAL_SPLIT);
	private final JSplitPane spectrumSplit=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	private final JSplitPane split=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	private final JTable table;
	private final TableRowSorter<TableModel> rowSorter;
	private final JTextField jtfFilter;
	private final DIAScanTableModel model;
	private final SearchParameters parameters;
	private final JTabbedPane primaryTabs=new JTabbedPane();
	
	
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
				final JFrame f=new JFrame("DIA Browser");
				f.addWindowListener(new WindowAdapter() {
					public void windowClosing(WindowEvent e) {
						System.exit(0);
					}
				});

				SearchParameters params=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), DataAcquisitionType.DIA, false, true, false);
				f.getContentPane().add(new DIABrowserPanel(params), BorderLayout.CENTER);

				f.pack();
				f.setSize(new Dimension(1900, 1030)); // for 1920x1080
				f.setVisible(true);
			}
		});

		Logger.logLine("Launching DIA Browser");
	}

	public DIABrowserPanel(SearchParameters parameters) {
		super(new BorderLayout());
		this.parameters=parameters;
		
		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(new LabeledComponent("<p style=\"font-size:12px; font-family: Helvetica, sans-serif\"><b>Parameters", new JLabel()));
		
		rawFileChooser=new FileChooserPanel(null, "RAW File", StripeFileGenerator.getFilenameFilter(), true) {
			private static final long serialVersionUID=1L;

			@Override
			public void update(File... filename) {
				super.update(filename);
				if (filename!=null&&filename.length>0&&filename[0]!=null) {
					updateRaw(filename[0]);
				}
			}
		};
		options.add(rawFileChooser);
		
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
		
		JPanel left=new JPanel(new BorderLayout());
		left.add(options, BorderLayout.NORTH);
		left.add(new JScrollPane(table), BorderLayout.CENTER);

		setLayout(new BorderLayout());
		left.add(searchPanel, BorderLayout.SOUTH);

		primaryTabs.addTab("Scans", rawSplit);
		rawSplit.setBottomComponent(spectrumSplit);
		primaryTabs.addTab(INTENSITY_DISTRIBUTION_TITLE, distributionSplit);
		primaryTabs.addTab(BOXPLOT_TITLE, boxplotSplit);
        
		split.setLeftComponent(left);
		split.setRightComponent(primaryTabs);
		
		table.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			@Override
			public void valueChanged(ListSelectionEvent e) {
				updateToSelected();
			}
		});
		
		add(split, BorderLayout.CENTER);
	}
	public void askForRaw() {
		rawFileChooser.askForFiles();
	}
	
	private float maxTIC=0.0f;
	private XYTrace chromatogram=null;
	private XYTrace precursorIntensityHistogram=null;
	private XYTrace fragmentIntensityHistogram=null;

	public void updateRaw(final File f) {
		SwingWorkerProgress<ArrayList<AcquiredSpectrum>> worker=new SwingWorkerProgress<ArrayList<AcquiredSpectrum>>((Frame)SwingUtilities.getWindowAncestor(this), "Please wait...", "Reading Raw File") {
			@Override
			protected ArrayList<AcquiredSpectrum> doInBackgroundForReal() throws Exception {
				ChartPanel structureChart=MzmlStructureCharter.getStructureChart(f);
				primaryTabs.addTab(STRUCTURE_TITLE, structureChart);
				primaryTabs.setSelectedIndex(primaryTabs.getTabCount()-1);
				
				dia=StripeFileGenerator.getFile(f, parameters);
				Logger.logLine("Read "+dia.getOriginalFileName()+", ("+dia.getRanges().size()+" total windows)");
				ArrayList<AcquiredSpectrum> scans=new ArrayList<AcquiredSpectrum>();
				Collection<XYPoint> tics=new ArrayList<XYPoint>();
				Collection<XYPoint> basepeaks=new ArrayList<XYPoint>();
				maxTIC=0.0f;
				
				ArrayList<PrecursorScan> precursors=dia.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE);
				
				int increment=Math.max(1, precursors.size()/1000);
				int scanCount=0;
				float tic=0.0f;
				TFloatFloatHashMap precursorIonDistribution=new TFloatFloatHashMap();
				double minMZ=Double.MAX_VALUE;
				double maxMZ=-Double.MAX_VALUE;
				
				for (PrecursorScan precursorScan : precursors) {
					if (precursorScan.getIsolationWindowLower()<minMZ) minMZ=precursorScan.getIsolationWindowLower();
					if (precursorScan.getIsolationWindowUpper()>maxMZ) maxMZ=precursorScan.getIsolationWindowUpper();
					
					scans.add(precursorScan);
					tic+=precursorScan.getTIC();
					basepeaks.add(new XYPoint(precursorScan.getScanStartTime()/60f, General.max(precursorScan.getIntensityArray())));

					if (scanCount%increment==0) {
						tics.add(new XYPoint(precursorScan.getScanStartTime()/60f, tic));
						if (tic>maxTIC) {
							maxTIC=tic;
						}
						tic=0;
					}
					scanCount++;
					
					for (float intensity : precursorScan.getIntensityArray()) {
						float bin=((int)(10.0f*Log.protectedLog10(intensity)))/10.0f;
						precursorIonDistribution.adjustOrPutValue(bin, 1.0f, 1.0f);
					}
				}
				
				float threshold=maxTIC/20f;
				float minRT=Float.MAX_VALUE;
				float maxRT=0.0f;
				for (XYPoint point : tics) {
					if (point.y>threshold) {
						if (point.x<minRT) minRT=(float)point.x;
						if (point.x>maxRT) maxRT=(float)point.x;
					}
				}
				XYTrace basepeak=new XYTrace(basepeaks, GraphType.area, "Precursor Basepeak", new Color(255, 0, 0, 50), 2.0f);
				chromatogram=new XYTrace(tics, GraphType.area, "Precursor TIC");
				precursorIntensityHistogram=new XYTrace(precursorIonDistribution, GraphType.area, "Log10 Precursor Intensity Distribution");
				
				ArrayList<PolymerIon> polymerList=PolymerIon.getAllPolymerProducts(new Range(minMZ, maxMZ));
				PolymerIon[] polymerIons=polymerList.toArray(new PolymerIon[polymerList.size()]);
				HashMap<PolymerIon, XYTrace> polymerMap=ChromatogramExtractor.extractFragmentChromatograms(parameters.getPrecursorTolerance(), polymerIons, scans, null, GraphType.boldline);
				
				TDoubleDoubleHashMap fragmentIonDistribution=new TDoubleDoubleHashMap();

				@SuppressWarnings("rawtypes")
				HashMap<Comparable, TFloatArrayList> maxIITByRange=new HashMap<>();
				@SuppressWarnings("rawtypes")
				HashMap<Comparable, TFloatArrayList> maxIITByRT=new HashMap<>();
				for (FragmentScan stripe : dia.getStripes(new Range(-Float.MAX_VALUE, Float.MAX_VALUE), -Float.MAX_VALUE, Float.MAX_VALUE, false)) {
					scans.add(stripe);
					
					float rtInMin=stripe.getScanStartTime()/60f;
					if (rtInMin>minRT&&rtInMin<maxRT) {
						@SuppressWarnings("rawtypes")
						Comparable key=stripe.getRange();
						TFloatArrayList iits=maxIITByRange.get(key);
						if (iits==null) {
							iits=new TFloatArrayList();
							maxIITByRange.put(key, iits);
						}
	
						key=5f*Math.round(stripe.getScanStartTime()/300f);
						TFloatArrayList rts=maxIITByRT.get(key);
						if (rts==null) {
							rts=new TFloatArrayList();
							maxIITByRT.put(key, rts);
						}
						
						iits.add(stripe.getIonInjectionTime()*1000f);
						rts.add(stripe.getIonInjectionTime()*1000f);
					}
					
					stripe.getScanStartTime();
					for (float intensity : stripe.getIntensityArray()) {
						double bin=((int)(10.0*Log.protectedLog10(intensity)))/10.0;
						fragmentIonDistribution.adjustOrPutValue(bin, 1.0, 1.0);
					}
				}
				fragmentIntensityHistogram=new XYTrace(fragmentIonDistribution, GraphType.area, "Log10 Fragment Intensity Distribution");
				
				double maxBPThreshold=basepeak.getMaxY()*0.001;
				
				ArrayList<XYTrace> polymerTraceList=new ArrayList<XYTrace>();
				polymerTraceList.add(0, XYTrace.round(basepeak, 3/60f));
				for (XYTrace trace : polymerMap.values()) {
					XYPoint maxXY = trace.getMaxXY();
					if (maxXY.y>maxBPThreshold) {
						polymerTraceList.add(XYTrace.round(trace, 3/60f));
						ArrayList<XYPoint> point=new ArrayList<>();
						point.add(maxXY);
						polymerTraceList.add(new XYTrace(point, GraphType.text, trace.getName(), trace.getColor(), Optional.ofNullable(null)));
					}
				}
				XYTrace[] array = polymerTraceList.toArray(new XYTrace[polymerTraceList.size()]);
				final ChartPanel contaminantIntensities=Charter.getChart("Log10 Contaminant Intensity", "Count", false, array);
				
				JTabbedPane tabs=new JTabbedPane();
				final ChartPanel precursorIntensities=Charter.getChart("Log10 Precursor Intensity", "Count", false, precursorIntensityHistogram);
				final ChartPanel fragmentIntensities=Charter.getChart("Log10 Fragment Intensity", "Count", false, fragmentIntensityHistogram);
				tabs.addTab("Precursor Intensity", precursorIntensities);
				tabs.addTab("Fragment Intensity", fragmentIntensities);
				distributionSplit.setTopComponent(contaminantIntensities);
				distributionSplit.setBottomComponent(tabs);
				distributionSplit.setDividerLocation(400);
				
				final ChartPanel iits=Charter.getBoxplotChart(null, "Precursor Isolation Window", "Ion Injection Time (in msec)", maxIITByRange);
				final ChartPanel rts=Charter.getBoxplotChart(null, "Retention Time Bin (in min)", "Ion Injection Time (in msec)", maxIITByRT);
				boxplotSplit.setTopComponent(iits);
				boxplotSplit.setBottomComponent(rts);
				boxplotSplit.setDividerLocation(400);
				
				Collections.sort(scans, new SpectrumComparator(SpectrumComparator.compareWithRT));
				
				return scans;
			}
			@Override
			protected void doneForReal(ArrayList<AcquiredSpectrum> t) {
				model.updateEntries(t);
				table.addRowSelectionInterval(0, 0);
			}
		};
		worker.execute();
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

		primaryTabs.setSelectedIndex(0);
	}

	public void resetScan(ArrayList<AcquiredSpectrum> entries) {
		int location=split.getDividerLocation();
		//System.out.println("location:"+location);
		if (location<=5) {
			location=400;
		}
		int locationRaw=rawSplit.getDividerLocation();
		//System.out.println("locationRaw:"+locationRaw);
		if (locationRaw<=5) {
			locationRaw=400;
		}
		int locationSpectrum=spectrumSplit.getDividerLocation();
		//System.out.println("locationRaw:"+locationRaw);
		if (locationSpectrum<=5) {
			locationSpectrum=400;
		}
		
		if (entries==null) {
			if (chromatogram!=null) {
				rawSplit.setTopComponent(Charter.getChart("Retention Time", "Precursor TIC", false, chromatogram));
			} else {
				split.setLeftComponent(new JLabel("Select a scan!"));
				return;
			}
		} else {
			final Spectrum spectrum;
			final double[] rtRange;
			if (entries.size()==1) {
				spectrum=entries.get(0);
				float rt=spectrum.getScanStartTime()/60f;
				rtRange=new double[] {rt, rt};
			} else {
				spectrum=SpectrumUtils.mergeSpectra(downcast(entries), parameters.getFragmentTolerance());
				float minRT=Float.MAX_VALUE;
				float maxRT=-Float.MAX_VALUE;
				for (Spectrum entry : entries) {
					float rt=entry.getScanStartTime()/60f;
					if (rt>maxRT) maxRT=rt;
					if (rt<minRT) minRT=rt;
				}
				rtRange=new double[] {maxRT, minRT};
			}

			final ChartPanel spectrumChart=Charter.getChart(spectrum);
			XYTrace intensityHistogram=new XYTrace(PivotTableGenerator.createPivotTable(Log.log10(spectrum.getIntensityArray())), GraphType.area, "Log10 Fragment Intensity Distribution");
			final ChartPanel precursorIntensities=Charter.getChart("Log10 Intensity", "Count (N="+spectrum.getIntensityArray().length+")", false, intensityHistogram);
			
			spectrumSplit.setLeftComponent(spectrumChart);
			spectrumSplit.setRightComponent(precursorIntensities);
			
			XYTrace marker=new XYTrace(rtRange, new double[] {0, maxTIC}, GraphType.dashedline, "marker");
			rawSplit.setTopComponent(Charter.getChart("Retention Time", "Precursor TIC", false, chromatogram, marker));
		}
		spectrumSplit.setDividerLocation(locationSpectrum);
		rawSplit.setDividerLocation(locationRaw);
		split.setDividerLocation(location);
	}
	
	private ArrayList<Spectrum> downcast(ArrayList<AcquiredSpectrum> spectra) {
		ArrayList<Spectrum> ret=new ArrayList<>();
		ret.addAll(spectra);
		return ret;
	}
}
