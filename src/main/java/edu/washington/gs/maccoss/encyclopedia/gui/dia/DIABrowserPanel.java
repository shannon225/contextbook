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

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
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

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LabeledComponent;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingWorkerProgress;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SpectrumComparator;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SpectrumUtils;

public class DIABrowserPanel extends JPanel {
	private static final long serialVersionUID=1L;
	public static final Color[] colors=new Color[] {Color.red, Color.blue, Color.green, Color.cyan, Color.magenta, Color.orange, Color.yellow, Color.pink, Color.gray, 
			Color.red.darker(), Color.blue.darker(), Color.green.darker(), Color.cyan.darker(), Color.magenta.darker(), Color.orange.darker(), Color.yellow.darker(), Color.pink.darker(), Color.gray.darker()};

	private final FileChooserPanel rawFileChooser;
	private final JSplitPane rawSplit=new JSplitPane(JSplitPane.VERTICAL_SPLIT);
	private final JSplitPane split=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	private final JTable table;
	private final TableRowSorter<TableModel> rowSorter;
	private final JTextField jtfFilter;
	private final DIAScanTableModel model;
	private final SearchParameters parameters;
	
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

				SearchParameters params=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), DataAcquisitionType.OVERLAPPING_DIA);
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
		
		rawFileChooser=new FileChooserPanel(null, "Raw", new SimpleFilenameFilter(".dia", ".mzML"), true) {
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

        
		split.setLeftComponent(left);
		split.setRightComponent(rawSplit);
		
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

	public void updateRaw(final File f) {
		SwingWorkerProgress<ArrayList<Spectrum>> worker=new SwingWorkerProgress<ArrayList<Spectrum>>((Frame)SwingUtilities.getWindowAncestor(this), "Please wait...", "Reading Raw File") {
			@Override
			protected ArrayList<Spectrum> doInBackgroundForReal() throws Exception {
				dia=MzmlToDIAConverter.getFile(f, parameters);
				Logger.logLine("Read "+dia.getOriginalFileName()+", ("+dia.getRanges().size()+" total windows)");
				ArrayList<Spectrum> scans=new ArrayList<Spectrum>();
				Collection<XYPoint> tics=new ArrayList<XYPoint>();
				maxTIC=0.0f;
				
				ArrayList<PrecursorScan> precursors=dia.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE);
				int increment=precursors.size()/1000;
				int scanCount=0;
				float tic=0.0f;
				for (PrecursorScan precursorScan : precursors) {
					scans.add(precursorScan);
					tic+=precursorScan.getTIC();
					if (scanCount%increment==0) {
						tics.add(new XYPoint(precursorScan.getScanStartTime()/60f, tic));
						if (tic>maxTIC) {
							maxTIC=tic;
						}
						tic=0;
					}
					scanCount++;
				}
				chromatogram=new XYTrace(tics, GraphType.area, "Precursor TIC");
				
				for (Stripe stripe : dia.getStripes(new Range(-Float.MAX_VALUE, Float.MAX_VALUE), -Float.MAX_VALUE, Float.MAX_VALUE, false)) {
					scans.add(stripe);
				}
				
				Collections.sort(scans, new SpectrumComparator(SpectrumComparator.compareWithRT));
				
				return scans;
			}
			@Override
			protected void doneForReal(ArrayList<Spectrum> t) {
				model.updateEntries(t);
			}
		};
		worker.execute();
	}

	public void updateToSelected() {
		int[] selection=table.getSelectedRows();
		if (selection.length<=0) return;
		
		ArrayList<Spectrum> entries=new ArrayList<Spectrum>();
		for (int row : selection) {
			Spectrum entry=model.getSelectedRow(table.convertRowIndexToModel(row));
			entries.add(entry);
		}
		resetScan(entries);
	}

	public void resetScan(ArrayList<Spectrum> entries) {
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
		
		
		if (entries==null) {
			if (chromatogram!=null) {
				rawSplit.setTopComponent(Charter.getChart("Retention Time", "Precursor TIC", false, chromatogram));
			} else {
				split.setLeftComponent(new JLabel("Select a scan!"));
				return;
			}
		} else {
				if (entries.size()==1) {
					rawSplit.setBottomComponent(Charter.getChart(entries.get(0)));
					float rt=entries.get(0).getScanStartTime()/60f;
					XYTrace marker=new XYTrace(new double[] {rt, rt}, new double[] {0, maxTIC}, GraphType.dashedline, "marker");
					rawSplit.setTopComponent(Charter.getChart("Retention Time", "Precursor TIC", false, chromatogram, marker));
				} else {
					rawSplit.setBottomComponent(Charter.getChart(SpectrumUtils.mergeSpectra(entries, parameters.getFragmentTolerance())));
					
					float minRT=Float.MAX_VALUE;
					float maxRT=-Float.MAX_VALUE;
					for (Spectrum spectrum : entries) {
						float rt=spectrum.getScanStartTime()/60f;
						if (rt>maxRT) maxRT=rt;
						if (rt<minRT) minRT=rt;
					}
					XYTrace marker=new XYTrace(new double[] {maxRT, minRT}, new double[] {0, maxTIC}, GraphType.dashedline, "marker");
					rawSplit.setTopComponent(Charter.getChart("Retention Time", "Precursor TIC", false, chromatogram, marker));
				}
		}
		rawSplit.setDividerLocation(locationRaw);
		split.setDividerLocation(location);
	}
	
	
}
