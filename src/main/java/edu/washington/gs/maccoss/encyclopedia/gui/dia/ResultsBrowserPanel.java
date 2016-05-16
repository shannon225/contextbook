package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.awt.BorderLayout;
import java.awt.Frame;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Optional;

import javax.swing.BoxLayout;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.SwingUtilities;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.algorithms.FragmentationScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.FragmentationTraceTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.LibraryPredictedFragmentationScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideQuantExtractorTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneFragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanRawScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LabeledComponent;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingWorkerProgress;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;

public class ResultsBrowserPanel extends JPanel {
	private static final long serialVersionUID=1L;

	private final FileChooserPanel blibFileChooser;
	private final FileChooserPanel rawFileChooser;
	private final JSplitPane dataSplit=new JSplitPane(JSplitPane.VERTICAL_SPLIT);
	private final JSplitPane rawSplit=new JSplitPane(JSplitPane.VERTICAL_SPLIT);
	private final JSplitPane peakPickingSplit=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	private final JSplitPane split=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	private final JTable table;
	private final LibraryEntryTableModel model;
	private final SearchParameters parameters;
	private final PecanRawScorer scorer;
	
	private StripeFileInterface dia=null;

	public ResultsBrowserPanel() {
		super(new BorderLayout());
		
		HashMap<String, String> map=SearchParameterParser.getDefaultParameters();
		map.put("-runPhosphoLocalization", "true");
		map.put("-deconvoluteOverlappingWindows", "true");
		map.put("-fixed", "");
		parameters=SearchParameterParser.parseParameters(map);
		scorer=new PecanRawScorer(parameters.getFragmentTolerance(), new LibraryPredictedFragmentationScorer(parameters));
		
		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(new LabeledComponent("<p style=\"font-size:12px; font-family: Helvetica, sans-serif\"><b>Parameters", new JLabel()));
		blibFileChooser=new FileChooserPanel(null, "Library", new SimpleFilenameFilter(".elib", ".blib"), true) {
			private static final long serialVersionUID=1L;

			@Override
			public void update(File... filename) {
				super.update(filename);
				if (filename!=null&&filename.length>0&&filename[0]!=null) {
					updateTable(filename[0]);
				}
			}
		};
		options.add(blibFileChooser);
		
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
		
		model=new LibraryEntryTableModel();
		table=new JTable(model) {
			private static final long serialVersionUID=1L;

			@Override
			public Object getValueAt(int row, int column) {
				if (column==0) return row+1;
				return super.getValueAt(row, column);
			}
		};
		table.setAutoCreateRowSorter(true);
		
		JPanel left=new JPanel(new BorderLayout());
		left.add(options, BorderLayout.NORTH);
		left.add(new JScrollPane(table), BorderLayout.CENTER);
		split.setLeftComponent(left);
		split.setRightComponent(dataSplit);
		
		table.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			@Override
			public void valueChanged(ListSelectionEvent e) {
				updateToSelected();
			}
		});
		
		add(split, BorderLayout.CENTER);
	}
	
	public void updateTable(final File f) {
		SwingWorkerProgress<ArrayList<LibraryEntry>> worker=new SwingWorkerProgress<ArrayList<LibraryEntry>>((Frame)SwingUtilities.getWindowAncestor(this), "Please wait...", "Reading Library") {
			@Override
			protected ArrayList<LibraryEntry> doInBackgroundForReal() throws Exception {
				LibraryInterface library=BlibToLibraryConverter.getFile(f);
				ArrayList<LibraryEntry> entries=library.getEntries(new Range(-Float.MAX_VALUE, Float.MAX_VALUE), false);
				return entries;
			}
			@Override
			protected void doneForReal(ArrayList<LibraryEntry> t) {
				model.updateEntries(t);
			}
		};
		worker.execute();
	}

	public void updateRaw(final File f) {
		SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame)SwingUtilities.getWindowAncestor(this), "Please wait...", "Reading Raw File") {
			@Override
			protected Nothing doInBackgroundForReal() throws Exception {
				dia=MzmlToDIAConverter.getFile(f, parameters);
				Logger.logLine("Read "+dia.getOriginalFileName()+", ("+dia.getRanges().size()+" total windows)");
				return Nothing.NOTHING;
			}
			@Override
			protected void doneForReal(Nothing t) {
				updateToSelected();
			}
		};
		worker.execute();
	}

	public void updateToSelected() {
		int[] selection=table.getSelectedRows();
		if (selection.length<=0) return;
		
		LibraryEntry entry=model.getSelectedRow(table.convertRowIndexToModel(selection[0]));
		resetPeptide(entry);
	}

	public void resetPeptide(LibraryEntry entry) {
		int location=dataSplit.getDividerLocation();
		System.out.println("location:"+location);
		if (location<=5) {
			location=800;
		}
		int locationRaw=rawSplit.getDividerLocation();
		System.out.println("locationRaw:"+locationRaw);
		if (locationRaw<=5) {
			locationRaw=400;
		}
		int locationPP=peakPickingSplit.getDividerLocation();
		System.out.println("locationPP:"+locationPP);
		if (locationPP<=5) {
			locationPP=400;
		}
		
		if (entry==null) {
			dataSplit.setTopComponent(new JLabel("Select a peptide!"));
			dataSplit.setBottomComponent(new JPanel());
			return;
		} else if (dia==null) {
			dataSplit.setTopComponent(new JLabel("Select a raw file!"));
			dataSplit.setBottomComponent(Charter.getChart(entry));
		} else {
			Logger.logLine("Parsing peptide...");
			PecanOneFragmentationModel model=new PecanOneFragmentationModel(new FastaPeptideEntry(entry.getPeptideModSeq()), parameters.getAAConstants());
			ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
			LibraryEntry unit=model.getUnitSpectrum(dia.getOriginalFileName(), entry.getAccessions(), (byte)entry.getPrecursorCharge(), entry.getRetentionTime(), parameters, 200.0);
			entries.add(unit);
			
			try {
				float rtRange=300f;
				ArrayList<Stripe> stripes=dia.getStripes(entry.getPrecursorMZ(), entry.getRetentionTime()-rtRange, entry.getRetentionTime()+rtRange, false);
				FragmentationTraceTask task=new FragmentationTraceTask(scorer, FragmentationTraceTask.PLOT_INTENSITIES, entries, stripes, new PrecursorScanMap(new ArrayList<PrecursorScan>()), parameters.getAAConstants());
				HashMap<LibraryEntry, PeptideScoringResult> result=task.call();
				
				ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
				for (Entry<LibraryEntry, PeptideScoringResult> resultEntry : result.entrySet()) {
					FragmentationScoringResult peptideResult=(FragmentationScoringResult)resultEntry.getValue();

					for (XYTrace trace : peptideResult.getFragmentationTraces()) {
						traces.add(trace);
					}
				}
				ChartPanel chart=Charter.getChart("RT ("+entry.getPrecursorMZ()+" M/Z)", "Intensity", true, traces.toArray(new XYTrace[traces.size()]));
				rawSplit.setTopComponent(chart);
				
				PSMData psmdata=new PSMData(entry.getAccessions(), entry.getSpectrumIndex(), entry.getPrecursorMZ(), entry.getPrecursorCharge(), entry.getPeptideModSeq(), entry.getRetentionTime(), entry.getScore(), 1.0f-entry.getScore(), 2*rtRange);
				PeptideQuantExtractorTask quantTask=new PeptideQuantExtractorTask(dia.getOriginalFileName(), psmdata, Optional.ofNullable((LibraryInterface)null), stripes, parameters, false);
				TransitionRefinementData data=quantTask.extractSpectrum(unit, 2*rtRange, false);
				if (parameters.isRunPhosphoLocalization()) {
					quantTask.runLocalization();
				}
				if (data!=null) {
					HashMap<String, ChartPanel> panels=TransitionRefiner.getChartPanels(data);
					peakPickingSplit.setLeftComponent(panels.get("median"));
					peakPickingSplit.setRightComponent(panels.get("unnormalized"));
					rawSplit.setBottomComponent(peakPickingSplit);
					peakPickingSplit.setDividerLocation(locationPP);
				} else {
					rawSplit.setBottomComponent(new JLabel("No quant data?"));
				}

				rawSplit.setDividerLocation(locationRaw);
				dataSplit.setTopComponent(rawSplit);
				
			} catch (Exception e) {
				JOptionPane.showMessageDialog(ResultsBrowserPanel.this, "Sorry, there was a problem reading the precursor window that contains ["+entry.getPrecursorMZ()+"]: "+e.getMessage(), "Error Reading DIA File",
						JOptionPane.ERROR_MESSAGE);
				e.printStackTrace();
			}
			
			dataSplit.setBottomComponent(Charter.getChart(entry));
			Logger.logLine("Finished reading peptide "+entry.getSpectrumName()+" (rt="+ entry.getRetentionTime()+")");
		}
		dataSplit.setDividerLocation(location);
	}
}
