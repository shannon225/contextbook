package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Frame;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.TreeMap;

import javax.swing.BoxLayout;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.SwingUtilities;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.axis.ValueAxis;

import edu.washington.gs.maccoss.encyclopedia.algorithms.FragmentationScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.FragmentationTraceTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.LibraryPredictedFragmentationScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideQuantExtractorTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneFragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanRawScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
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
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LabeledComponent;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingWorkerProgress;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class ResultsBrowserPanel extends JPanel {
	private static final long serialVersionUID=1L;
	public static final Color[] colors=new Color[] {Color.red, Color.blue, Color.green, Color.cyan, Color.magenta, Color.orange, Color.yellow, Color.pink, Color.gray};

	private final FileChooserPanel elibFileChooser;
	private final FileChooserPanel rawFileChooser;
	private final JSplitPane dataSplit=new JSplitPane(JSplitPane.VERTICAL_SPLIT);
	private final JSplitPane rawSplit=new JSplitPane(JSplitPane.VERTICAL_SPLIT);
	private final JSplitPane peakPickingSplit=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	private final JSplitPane split=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	private final JTable table;
	private final LibraryEntryTableModel model;
	private final SearchParameters parameters;
	private final PecanRawScorer scorer;
	
	private LibraryInterface library=null;
	private StripeFileInterface dia=null;

	public ResultsBrowserPanel(SearchParameters parameters) {
		super(new BorderLayout());
		this.parameters=parameters;
		scorer=new PecanRawScorer(parameters.getFragmentTolerance(), new LibraryPredictedFragmentationScorer(parameters));
		
		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(new LabeledComponent("<p style=\"font-size:12px; font-family: Helvetica, sans-serif\"><b>Parameters", new JLabel()));
		elibFileChooser=new FileChooserPanel(null, "Library", new SimpleFilenameFilter(".elib"), true) {
			private static final long serialVersionUID=1L;

			@Override
			public void update(File... filename) {
				super.update(filename);
				if (filename!=null&&filename.length>0&&filename[0]!=null) {
					updateTable(filename[0]);
				}
			}
		};
		options.add(elibFileChooser);
		
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
				library=BlibToLibraryConverter.getFile(f);
				ArrayList<LibraryEntry> entries=library.getEntries(new Range(-Float.MAX_VALUE, Float.MAX_VALUE), false);
				
				Optional<StripeFileInterface> source=library.getSource(parameters);
				if (source.isPresent()) {
					dia=source.get();
				}
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
				float rtRange=30f;
				ArrayList<Stripe> stripes=dia.getStripes(entry.getPrecursorMZ(), entry.getRetentionTime()-rtRange, entry.getRetentionTime()+rtRange, false);
				FragmentationTraceTask task=new FragmentationTraceTask(scorer, FragmentationTraceTask.PLOT_INTENSITIES, entries, stripes, new PrecursorScanMap(new ArrayList<PrecursorScan>()), parameters.getAAConstants());
				HashMap<LibraryEntry, PeptideScoringResult> result=task.call();
				
				ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
				for (Entry<LibraryEntry, PeptideScoringResult> resultEntry : result.entrySet()) {
					FragmentationScoringResult peptideResult=(FragmentationScoringResult)resultEntry.getValue();

					for (XYTrace trace : peptideResult.getFragmentationTraces()) {
						double[] intensities=trace.toArrays().y;
						for (int i=0; i<intensities.length; i++) {
							if (intensities[i]>0) {
								//traces.add(SkylineSGFilter.paddedSavitzkyGolaySmooth(trace));
								traces.add(trace);
								break;
							}
						}
					}
				}
				
				
				/*for (XYTrace xyTrace : traces) {
					System.out.println("float[] "+xyTrace.getName()+"=new float[] {");
					double[] d=xyTrace.toArrays().y;
					for (int i=0; i<d.length; i++) {
						System.out.print(d[i]+"f, ");
					}
					System.out.println("};");
					System.out.println("chromatograms.add("+xyTrace.getName()+");");
				}

				System.out.print("String[] ionNames=new String[] {");
				for (XYTrace xyTrace : traces) {
					System.out.print("\""+xyTrace.getName()+"\", ");
				}
				System.out.println("};");
				System.out.println("float[] rts=new float[] {");
				for (double d : traces.get(0).toArrays().x) {
					System.out.print(d+"f, ");
				}
				System.out.println("};");
				*/
				
				ChartPanel chart=Charter.getChart("Retention Time", "Intensity", true, traces.toArray(new XYTrace[traces.size()]));
				rawSplit.setTopComponent(chart);

				PhosphoLocalizer localizer=new PhosphoLocalizer(dia, library, parameters);
				PSMData psmdata=new PSMData(entry.getAccessions(), entry.getSpectrumIndex(), entry.getPrecursorMZ(), entry.getPrecursorCharge(), entry.getPeptideModSeq(), entry.getRetentionTime(), entry.getScore(), 1.0f-entry.getScore(), 2*rtRange);
				PeptideQuantExtractorTask quantTask=new PeptideQuantExtractorTask(dia.getOriginalFileName(), psmdata, Optional.ofNullable(localizer), stripes, parameters, false);
				TransitionRefinementData data=quantTask.extractSpectrum(unit, 2*rtRange, false);
				if (data!=null) {
					HashMap<String, ChartPanel> panels=TransitionRefiner.getChartPanels(data);
					peakPickingSplit.setLeftComponent(panels.get("median"));
					peakPickingSplit.setRightComponent(panels.get("unnormalized"));
					
					JTabbedPane tabs=new JTabbedPane();
					
					Optional<PhosphoLocalizationData> phosphoData=Optional.empty();
					if (parameters.isRunPhosphoLocalization()) {
						phosphoData=quantTask.runLocalization();
					}
					if (phosphoData.isPresent()) {
						HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> allVsUniqueList=phosphoData.get().getTraces();
						HashMap<String, XYTrace[]> uniqueFragmentIons=phosphoData.get().getUniqueFragmentIons();
						HashMap<String, XYPoint> localizationScores=phosphoData.get().getLocalizationScores();
						
						ArrayList<XYTrace> phosphoTraces=new ArrayList<XYTrace>();

						TreeMap<String, ChartPanel> panelMap=new TreeMap<String, ChartPanel>();
						int i=0;
						for (Entry<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> phosphoentry : allVsUniqueList.entrySet()) {
							String seq=phosphoentry.getKey();
							Pair<TFloatFloatHashMap, TFloatFloatHashMap> pair=phosphoentry.getValue();
							Color color=i>=colors.length?colors[i-colors.length].brighter():colors[i];
							//phosphoTraces.add(new XYTrace(pair.x, GraphType.line, "ALL_"+seq, new Color(color.getRed(), color.getGreen(), color.getBlue(), 150), 4.0f));
							phosphoTraces.add(new XYTrace(pair.y, GraphType.line, "UNI_"+seq, color, 2.0f));

							XYPoint point=localizationScores.get(seq);
							phosphoTraces.add(new XYTrace(new double[] {point.x}, new double[] {point.y}, GraphType.point, "center", color, 2.0f));
							i++;
							
							XYTrace[] uniqueFragments=uniqueFragmentIons.get(seq);
							
							panelMap.put(seq+" ("+(Math.round(point.y*10.0f)/10.0f)+")", Charter.getChart("Retention Time", "Intensity", true, uniqueFragments));
						}
						ChartPanel phosphoPane=Charter.getChart("Retention Time", "Score", true, phosphoTraces.toArray(new XYTrace[phosphoTraces.size()]));
						ValueAxis axis=phosphoPane.getChart().getXYPlot().getRangeAxis();
						org.jfree.data.Range range=axis.getRange();
						axis.setRange(new org.jfree.data.Range(0.0f, Math.max(2.0f, range.getUpperBound())));
						tabs.add("Phospho Localization", phosphoPane);
						tabs.add("Unique Fragment Ions", Charter.getTabbedChartPane(panelMap));
					}
					tabs.add("Quantification", peakPickingSplit);
					
					rawSplit.setBottomComponent(tabs);
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
