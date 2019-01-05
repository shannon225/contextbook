package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Frame;
import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Optional;

import javax.swing.BoxLayout;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
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
import org.jfree.chart.axis.ValueAxis;

import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.PeptideQuantExtractorTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SimplePeptidePrecursor;
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
import edu.washington.gs.maccoss.encyclopedia.gui.massspec.FragmentationTable;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.ChromatogramExtractor;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class ResultsBrowserPanel extends JPanel {
	private static final long serialVersionUID=1L;
	public static final Color[] colors=new Color[] {Color.red, Color.blue, Color.green, Color.cyan, Color.magenta, Color.orange, Color.yellow, Color.pink, Color.gray, 
			Color.red.darker(), Color.blue.darker(), Color.green.darker(), Color.cyan.darker(), Color.magenta.darker(), Color.orange.darker(), Color.yellow.darker(), Color.pink.darker(), Color.gray.darker()};

	private final FileChooserPanel elibFileChooser;
	private final FileChooserPanel rawFileChooser;
	private final JSplitPane dataSplit=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	private final JSplitPane rawSplit=new JSplitPane(JSplitPane.VERTICAL_SPLIT);
	private final JSplitPane peakPickingSplit=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	private final JSplitPane split=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	private final JTable table;
	private final TableRowSorter<TableModel> rowSorter;
	private final JTextField jtfFilter;
	private final LibraryEntryTableModel model;
	private final SearchParameters parameters;
	private final float minimumScore;
	
	private LibraryInterface library=null;
	private StripeFileInterface dia=null;
	private Optional<PhosphoLocalizer> nullableLocalizer=Optional.ofNullable(null);


	public ResultsBrowserPanel(SearchParameters parameters) {
		super(new BorderLayout());
		this.parameters=parameters;
		this.minimumScore=-Log.log10(parameters.getPercolatorThreshold());
		
		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(new LabeledComponent("<p style=\"font-size:12px; font-family: Helvetica, sans-serif\"><b>Parameters", new JLabel()));
		elibFileChooser=new FileChooserPanel(null, "Library", new SimpleFilenameFilter(LibraryFile.DLIB, LibraryFile.ELIB), true) {
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
		
		rawFileChooser=new FileChooserPanel(null, "Raw", StripeFileGenerator.getFilenameFilter(), true) {
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
		left.add(new JScrollPane(table, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED), BorderLayout.CENTER);

		setLayout(new BorderLayout());
		left.add(searchPanel, BorderLayout.SOUTH);
		left.setMinimumSize(new Dimension(100, 100));
        
		split.setLeftComponent(left);
		split.setRightComponent(new JLabel("Select a peptide!"));
		
		table.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			@Override
			public void valueChanged(ListSelectionEvent e) {
				updateToSelected();
			}
		});
		
		add(split, BorderLayout.CENTER);
	}
	
	public void askForLibrary() {
		elibFileChooser.askForFiles();
	}
	public void askForRaw() {
		rawFileChooser.askForFiles();
	}
	
	public void updateTable(final File f) {
		SwingWorkerProgress<ArrayList<LibraryEntry>> worker=new SwingWorkerProgress<ArrayList<LibraryEntry>>((Frame)SwingUtilities.getWindowAncestor(this), "Please wait...", "Reading Library") {
			@Override
			protected ArrayList<LibraryEntry> doInBackgroundForReal() throws Exception {
				LibraryFile.OPEN_IN_PLACE=true;
				library=BlibToLibraryConverter.getFile(f);
				LibraryFile.OPEN_IN_PLACE=false;

				ArrayList<LibraryEntry> entries=library.getEntries(new Range(-Float.MAX_VALUE, Float.MAX_VALUE), false, parameters.getAAConstants());

				final Optional<Path> source = library.getSource(parameters);
				if (source.isPresent()) {
					try {
						dia = StripeFileGenerator.getFile(source.get().toFile(), parameters); // assumes the .DIA file exists or should be created
					} catch (Exception e) {
						Logger.errorLine("Sorry, can't load DIA from library annotation. Please load it after!");
						Logger.errorException(e);
					}
				}

				if (dia!=null&&library!=null&&parameters.getLocalizingModification().isPresent()) {
					PhosphoLocalizer localizer=new PhosphoLocalizer(dia, parameters.getLocalizingModification().get(), library, parameters);
					nullableLocalizer=Optional.ofNullable(localizer);
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
				
				dia=StripeFileGenerator.getFile(f, parameters, true);

				if (dia!=null&&library!=null&&parameters.getLocalizingModification().isPresent()) {
					PhosphoLocalizer localizer=new PhosphoLocalizer(dia, parameters.getLocalizingModification().get(), library, parameters);
					nullableLocalizer=Optional.ofNullable(localizer);
				}

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
		int locationPP=peakPickingSplit.getDividerLocation();
		//System.out.println("locationPP:"+locationPP);
		if (locationPP<=5) {
			locationPP=400;
		}
		
		if (entry==null) {
			split.setLeftComponent(new JLabel("Select a peptide!"));
			return;
		} else if (dia==null) {
			dataSplit.setLeftComponent(Charter.getChart(new AnnotatedLibraryEntry(entry, parameters)));
			dataSplit.setRightComponent(new FragmentationTable(entry, entry.getPeptideModSeq(), parameters));
			split.setRightComponent(dataSplit);
		} else {
			Logger.logLine("Parsing peptide...");
			FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
			ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
			float targetRT=entry.getRetentionTime();
			AnnotatedLibraryEntry unit=model.getUnitSpectrum(dia.getOriginalFileName(), entry.getAccessions(), (byte)entry.getPrecursorCharge(), targetRT, parameters, 0.0, true);
			entries.add(unit);
			
			try {
				float rtRange=parameters.getLocalizingModification().isPresent()?dia.getGradientLength()/20.0f:(3f*parameters.getExpectedPeakWidth());
				
				ArrayList<FragmentScan> stripes=dia.getStripes(entry.getPrecursorMZ(), targetRT-rtRange, targetRT+rtRange, false);
				Collections.sort(stripes);

				Float targetRTFloat=targetRT;
				ArrayList<Spectrum> downcastedSpectra=FragmentScan.downcastStripeToSpectrum(stripes);
				HashMap<FragmentIon, XYTrace> fragmentTraceMap=ChromatogramExtractor.extractFragmentChromatograms(parameters.getFragmentTolerance(), model.getPrimaryIonObjects(parameters.getFragType(), (byte)entry.getPrecursorCharge(), true), downcastedSpectra, targetRTFloat, GraphType.line);
				ArrayList<XYTrace> traces=new ArrayList<XYTrace>(fragmentTraceMap.values());
				Collections.sort(traces);

				ArrayList<Spectrum> precursors=PrecursorScan.downcast(dia.getPrecursors(targetRT-rtRange, targetRT+rtRange));
				ChartPanel precursorChart=Charter.getChart("Retention Time", "Intensity", true, ChromatogramExtractor.extractPrecursorChromatograms(parameters.getPrecursorTolerance(), entry.getPrecursorMZ(), entry.getPrecursorCharge(), precursors));
				
				ChartPanel fragmentChart=Charter.getChart("Retention Time (min)", "Intensity", true, traces.toArray(new XYTrace[traces.size()]));
				JTabbedPane primaryTabs=new JTabbedPane();
				primaryTabs.add("Fragments", fragmentChart);
				primaryTabs.add("Precursors", precursorChart);
				rawSplit.setTopComponent(primaryTabs);
				
				PSMData psmdata=new PSMData(entry.getAccessions(), entry.getSpectrumIndex(), entry.getPrecursorMZ(), entry.getPrecursorCharge(), entry.getPeptideModSeq(), targetRT, entry.getScore(), 1.0f-entry.getScore(), 2*rtRange, false, parameters.getAAConstants());
				PeptideQuantExtractorTask quantTask=new PeptideQuantExtractorTask(dia.getOriginalFileName(), psmdata, Optional.empty(), nullableLocalizer, stripes, parameters, false);
				TransitionRefinementData data=quantTask.extractSpectrum(unit, rtRange, false, false, false);
				if (data!=null) {
					HashMap<String, ChartPanel> panels=TransitionRefiner.getChartPanels(data);
					peakPickingSplit.setLeftComponent(panels.get("median"));
					peakPickingSplit.setRightComponent(panels.get("unnormalized"));
					
					JTabbedPane tabs=new JTabbedPane();
					
					Optional<PhosphoLocalizationData> maybePhosphoData=Optional.empty();
					if (parameters.getLocalizingModification().isPresent()) {
						maybePhosphoData=quantTask.runLocalization(true);
					}
					if (maybePhosphoData.isPresent()) {
						PhosphoLocalizationData actuallyPhosphoData=maybePhosphoData.get();
						HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> ionsVsUniqueList=actuallyPhosphoData.getScoreTraces();
						HashMap<String, HashMap<FragmentIon, XYTrace>> uniqueFragmentIons=actuallyPhosphoData.getUniqueFragmentIons();
						HashMap<String, HashMap<FragmentIon, XYTrace>> otherFragmentIons=actuallyPhosphoData.getOtherFragmentIons();
						
						HashMap<String, XYPoint> localizationScores=actuallyPhosphoData.getLocalizationScores();

						ArrayList<XYTrace> complementaryIonsTraces=new ArrayList<XYTrace>();
						ArrayList<XYTrace> phosphoTraces=new ArrayList<XYTrace>();

						JTabbedPane tabPanel=new JTabbedPane();
						HashMap<String, String> keyVsName=new HashMap<String, String>();
						int colorIndex=0;
						for (Entry<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> phosphoentry : ionsVsUniqueList.entrySet()) {
							String sequenceKey=phosphoentry.getKey();
							String peptideModSeq=sequenceKey.replaceAll("\\(", "").replaceAll("\\)", "");
							Pair<TFloatFloatHashMap, TFloatFloatHashMap> pair=phosphoentry.getValue();
							Color color=colorIndex>=colors.length?colors[colorIndex-colors.length].brighter():colors[colorIndex];
							colorIndex++;
							
							complementaryIonsTraces.add(new XYTrace(pair.x, GraphType.line, sequenceKey, color, 2.0f));
							phosphoTraces.add(new XYTrace(pair.y, GraphType.line, sequenceKey, color, 2.0f));
							
							XYPoint point=localizationScores.get(sequenceKey);
							if (point!=null) {
								//complementaryIonsTraces.add(new XYTrace(new double[] {point.x/60f}, new double[] {0}, GraphType.point, "center", color, 2.0f)); //FIXME
								//phosphoTraces.add(new XYTrace(new double[] {point.x/60f}, new double[] {point.y}, GraphType.point, "center", color, 2.0f)); //FIXME
							} else {
								continue;
							}
							
							HashMap<FragmentIon, XYTrace> uniqueFragments=uniqueFragmentIons.get(sequenceKey);
							HashMap<FragmentIon, XYTrace> otherFragments=new HashMap<FragmentIon, XYTrace>(otherFragmentIons.get(sequenceKey));
							for (FragmentIon ion : new ArrayList<FragmentIon>(otherFragments.keySet())) {
								if (ion.index<=1) {
									//otherFragments.remove(ion);
								}
							}
							
							HashMap<FragmentIon, XYTrace> allFragments=new HashMap<FragmentIon, XYTrace>();
							allFragments.putAll(uniqueFragments);
							allFragments.putAll(otherFragments);
							
							ArrayList<XYTrace> uniqueFragmentsList=new ArrayList<XYTrace>(allFragments.values());
							double maxPoint=XYTrace.getMaxY(uniqueFragmentsList);
							
							//uniqueFragmentsList.add(new XYTrace(new double[] {point.x/60f, point.x/60f}, new double[] {0.0, maxPoint}, GraphType.dashedline, "center", Color.BLACK, 2.0f)); // FIXME
							XYTraceInterface[] fragmentTraces=uniqueFragmentsList.toArray(new XYTrace[uniqueFragmentsList.size()]);
							
							if (point.y>=minimumScore&&actuallyPhosphoData.getPassingForms().containsKey(sequenceKey)) {
								keyVsName.put(sequenceKey, sequenceKey+" ("+(Math.round(point.y*10.0f)/10.0f)+")");
							} else {
								keyVsName.put(sequenceKey, sequenceKey+" (not sig)");
							}

							//ChartPanel chartPane=Charter.getChart("Retention Time (min)", "Intensity", true, fragmentTraces); // FIXME
							ChartPanel chartPane=Charter.getChart("Retention Time (min)", "Intensity", false, fragmentTraces);
							
							//TransitionRefinementData quantData=actuallyPhosphoData.getPassingForms().get(sequenceKey);
							//AnnotatedLibraryEntry annotatedEntry;
							/*if (quantData!=null) {
								annotatedEntry=quantData.getEntry(unit, parameters);
							} else {
								Spectrum bestStripe=ChromatogramExtractor.getTargetStripeByRT(downcastedSpectra, (float)point.x);
								annotatedEntry=new AnnotatedLibraryEntry(new SimplePeptidePrecursor(peptideModSeq, entry.getPrecursorCharge()), bestStripe, parameters);
							}*/
							Spectrum bestStripe=ChromatogramExtractor.getTargetStripeByRT(downcastedSpectra, (float)point.x);
							AnnotatedLibraryEntry annotatedEntry=new AnnotatedLibraryEntry(new SimplePeptidePrecursor(peptideModSeq, entry.getPrecursorCharge(), parameters.getAAConstants()), bestStripe, parameters);

							JPanel specFragPane=new JPanel(new BorderLayout());
							ChartPanel spectrumPane=Charter.getChart(annotatedEntry);
							FragmentationTable fragTable=new FragmentationTable(actuallyPhosphoData, sequenceKey, parameters);
							specFragPane.add(spectrumPane, BorderLayout.NORTH);
							specFragPane.add(fragTable, BorderLayout.CENTER);
							
							JPanel localizationPane=new JPanel(new BorderLayout());
							localizationPane.add(chartPane, BorderLayout.CENTER);
							localizationPane.add(specFragPane, BorderLayout.EAST);
							
							tabPanel.addTab(keyVsName.get(sequenceKey), localizationPane);
						}
						
						ChartPanel phosphoPane=Charter.getChart("Retention Time (min)", "Score", true, phosphoTraces.toArray(new XYTrace[phosphoTraces.size()]));
						ValueAxis axis=phosphoPane.getChart().getXYPlot().getRangeAxis();
						org.jfree.data.Range range=axis.getRange();
						axis.setRange(new org.jfree.data.Range(0.0f, Math.max(2.0f, range.getUpperBound())));
						tabs.add("Phospho Localization", phosphoPane);
						
						ChartPanel coelutingPane=Charter.getChart("Retention Time (min)", "Coeluting Ions", true, complementaryIonsTraces.toArray(new XYTrace[complementaryIonsTraces.size()]));
						axis=coelutingPane.getChart().getXYPlot().getRangeAxis();
						range=axis.getRange();
						axis.setRange(new org.jfree.data.Range(0.0f, Math.max(2.0f, range.getUpperBound())));
						tabs.add("Coeluting Ions", coelutingPane);
						
						tabs.add("Unique Fragment Ions", tabPanel);
						
					}
					tabs.add("Quantification", peakPickingSplit);
					if (!maybePhosphoData.isPresent()) {
						Spectrum bestStripe=ChromatogramExtractor.getTargetStripeByRT(downcastedSpectra, targetRT);
						AnnotatedLibraryEntry annotatedEntry=new AnnotatedLibraryEntry(new SimplePeptidePrecursor(entry.getPeptideModSeq(), entry.getPrecursorCharge(), parameters.getAAConstants()), bestStripe, parameters);

						JSplitPane specFragPane=new JSplitPane(JSplitPane.VERTICAL_SPLIT);
						ChartPanel spectrumPane=Charter.getChart(annotatedEntry);
						FragmentationTable fragTable=new FragmentationTable(entry, entry.getPeptideModSeq(), parameters);
						specFragPane.add(spectrumPane, JSplitPane.TOP);
						specFragPane.add(fragTable, JSplitPane.BOTTOM);

						tabs.add("Detection", specFragPane);
					}
					
					rawSplit.setBottomComponent(tabs);
					peakPickingSplit.setDividerLocation(locationPP);
				} else {
					rawSplit.setBottomComponent(new JLabel("No quant data?"));
				}

				rawSplit.setDividerLocation(locationRaw);
				split.setRightComponent(rawSplit);
				
			} catch (Exception e) {
				JOptionPane.showMessageDialog(ResultsBrowserPanel.this, "Sorry, there was a problem reading the precursor window that contains ["+entry.getPrecursorMZ()+"]: "+e.getMessage(), "Error Reading DIA File",
						JOptionPane.ERROR_MESSAGE);
				e.printStackTrace();
			}
			Logger.logLine("Finished reading peptide "+entry.getSpectrumName()+" (rt="+ targetRT+")");
		}
		split.setDividerLocation(location);
	}
}
