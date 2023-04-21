package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Frame;
import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
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
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LocalizedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
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
import edu.washington.gs.maccoss.encyclopedia.gui.massspec.FragmentationTable;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.AcquiredSpectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.ChromatogramExtractor;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;

public class LocalizationResultsBrowserPanel extends JPanel {
	private static final long serialVersionUID=1L;
	public static final Color[] colors=new Color[] {Color.red, Color.blue, Color.green, Color.cyan, Color.magenta, Color.orange, Color.yellow, Color.pink, Color.gray, 
			Color.red.darker(), Color.blue.darker(), Color.green.darker(), Color.cyan.darker(), Color.magenta.darker(), Color.orange.darker(), Color.yellow.darker(), Color.pink.darker(), Color.gray.darker()};

	private final FileChooserPanel elibFileChooser;
	private final FileChooserPanel rawFileChooser;
	private final JSplitPane tableDataSplit=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	private final JSplitPane chartSplit=new JSplitPane(JSplitPane.VERTICAL_SPLIT);
	private final JTable table;
	private final TableRowSorter<TableModel> rowSorter;
	private final JTextField jtfFilter;
	private final LocalizedLibraryEntryTableModel model;
	private final SearchParameters parameters;
	
	private LibraryInterface library=null;
	private StripeFileInterface dia=null;


	public LocalizationResultsBrowserPanel(SearchParameters parameters) {
		super(new BorderLayout());
		this.parameters=parameters;
		
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
		
		model=new LocalizedLibraryEntryTableModel();
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
		table.getColumnModel().getColumn(LocalizedLibraryEntryTableModel.deltaRTColumnIndex).setCellRenderer(new StatusColumnCellRenderer());

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
        
		tableDataSplit.setLeftComponent(left);
		tableDataSplit.setRightComponent(new JLabel("Select a peptide!"));
		
		table.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			@Override
			public void valueChanged(ListSelectionEvent e) {
				updateToSelected();
			}
		});
		
		add(tableDataSplit, BorderLayout.CENTER);
	}
	
	public void askForLibrary() {
		elibFileChooser.askForFiles();
	}
	public void askForRaw() {
		rawFileChooser.askForFiles();
	}
	
	public void updateTable(final File f) {
		SwingWorkerProgress<ArrayList<LocalizedLibraryEntry>> worker=new SwingWorkerProgress<ArrayList<LocalizedLibraryEntry>>((Frame)SwingUtilities.getWindowAncestor(this), "Please wait...", "Reading Library") {
			@Override
			protected ArrayList<LocalizedLibraryEntry> doInBackgroundForReal() throws Exception {
				ArrayList<LocalizedLibraryEntry> entries=new ArrayList<>();
				library=BlibToLibraryConverter.getFile(f);
				if (library instanceof LibraryFile) {
					entries=((LibraryFile)library).getAllLocalizedEntries(-100000f, false, parameters.getLocalizingModification().get(), false, parameters.getAAConstants());
				}
				
				if (entries.size()==0) {
					JOptionPane.showMessageDialog(LocalizationResultsBrowserPanel.this, "Sorry, you can only view localized results from Thesaurus", "Invalid results data", JOptionPane.WARNING_MESSAGE);
					return entries;
				}
				
				final Optional<Path> source = library.getSource(parameters);
				if (source.isPresent()) {
					try {
						dia = StripeFileGenerator.getFile(source.get().toFile(), parameters); // assumes the .DIA file exists or should be created
					} catch (EncyclopediaException ee) {
						Logger.errorException(ee);
					}
				}
				return entries;				
			}
			@Override
			protected void doneForReal(ArrayList<LocalizedLibraryEntry> t) {
				model.updateEntries(t);
			}
		};
		worker.execute();
	}

	public void updateRaw(final File f) {
		SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame)SwingUtilities.getWindowAncestor(this), "Please wait...", "Reading Raw File") {
			@Override
			protected Nothing doInBackgroundForReal() throws Exception {
				dia=StripeFileGenerator.getFile(f, parameters);

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
		
		ArrayList<LocalizedLibraryEntry> entry=model.getSelectedRow(table.convertRowIndexToModel(selection[0]));
		resetPeptide(entry);
	}

	public void resetPeptide(ArrayList<LocalizedLibraryEntry> entries) {
		int location=tableDataSplit.getDividerLocation();
		//System.out.println("location:"+location);
		if (location<=5) {
			location=400;
		}
		int locationRaw=chartSplit.getDividerLocation();
		//System.out.println("locationRaw:"+locationRaw);
		if (locationRaw<=5) {
			locationRaw=200;
		}
		
		if (entries==null) {
			tableDataSplit.setRightComponent(new JLabel("Select a peptide!"));
			return;
		} else if (dia==null) {
			JTabbedPane tabs=new JTabbedPane();
			for (LocalizedLibraryEntry entry : entries) {
				JSplitPane sp=new JSplitPane(JSplitPane.VERTICAL_SPLIT);
				sp.setTopComponent(Charter.getChart(new AnnotatedLibraryEntry(entry, parameters)));
				sp.setBottomComponent(new FragmentationTable(entry, entry.getPeptideModSeq(), parameters));
				tabs.add(entry.getPeptideAnnotation().getNonDirectionalPeptideAnnotation()+" ("+(Math.round(entry.getLocalizationScore()*10.0f)/10.0f)+")", sp);
			}
			tableDataSplit.setRightComponent(tabs);
		} else {
			float deltaRT=60f;
			float minRT=Float.MAX_VALUE;
			float maxRT=-Float.MAX_VALUE;
			for (LocalizedLibraryEntry entry : entries) {
				if (entry.getRetentionTime()>maxRT) maxRT=entry.getRetentionTime();
				if (entry.getRetentionTime()<minRT) minRT=entry.getRetentionTime();
			}
			double precursorMZ = entries.get(0).getPrecursorMZ();
			byte precursorCharge = entries.get(0).getPrecursorCharge();

			try {
				ArrayList<Spectrum> precursors=PrecursorScan.downcast(dia.getPrecursors(minRT-deltaRT, maxRT+deltaRT));
				ArrayList<Spectrum> trimmedPrecursors=new ArrayList<>();
				for (Spectrum spectrum : precursors) {
					if (spectrum instanceof AcquiredSpectrum) {
						AcquiredSpectrum acquiredSpectrum=(AcquiredSpectrum)spectrum;

						if (precursorMZ>acquiredSpectrum.getIsolationWindowLower()&&precursorMZ<acquiredSpectrum.getIsolationWindowUpper()) {
							trimmedPrecursors.add(acquiredSpectrum);
						}
					}
				}
				precursors=trimmedPrecursors;
				
				ArrayList<FragmentScan> stripes=dia.getStripes(precursorMZ, minRT-deltaRT, maxRT+deltaRT, false);
				
				XYTraceInterface[] precursorTraces = ChromatogramExtractor.extractPrecursorChromatograms(parameters.getPrecursorTolerance(), precursorMZ, precursorCharge, precursors, true, false);
				double maxPrecursor=XYTrace.getMaxY(precursorTraces);
				ArrayList<XYTraceInterface> precursorTraceList=new ArrayList<>();
				for (XYTraceInterface trace : precursorTraces) {
					precursorTraceList.add(trace);
				}
				for (LocalizedLibraryEntry entry : entries) {
					precursorTraceList.add(new XYTrace(new double[] {entry.getRetentionTime()/60f, entry.getRetentionTime()/60f}, new double[] {0.0,  maxPrecursor}, GraphType.dashedline, entry.getPeptideModSeq(), Color.darkGray, 2.0f));
				}
				
				ChartPanel precursorChart=Charter.getChart("Retention Time", "Precursor Intensity", false, precursorTraceList.toArray(new XYTraceInterface[precursorTraceList.size()]));
	
				tableDataSplit.setRightComponent(chartSplit);
				chartSplit.setTopComponent(precursorChart);
	
				JTabbedPane tabs=new JTabbedPane();
				for (LocalizedLibraryEntry entry : entries) {
					FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
					
					ArrayList<Spectrum> downcastedSpectra=new ArrayList<Spectrum>();
					for (Spectrum spectrum : stripes) {
						if (spectrum.getScanStartTime()>entry.getRetentionTime()-deltaRT&&spectrum.getScanStartTime()<entry.getRetentionTime()+deltaRT) {
							downcastedSpectra.add(spectrum);
						}
					}
					HashMap<FragmentIon, XYTrace> fragmentTraceMap=ChromatogramExtractor.extractFragmentChromatograms(parameters.getFragmentTolerance(), model.getPrimaryIonObjects(parameters.getFragType(), (byte)entry.getPrecursorCharge(), true, true), downcastedSpectra, entry.getRetentionTime(), GraphType.dashedline, true, false);
					HashMap<FragmentIon, XYTrace> targetMap=ChromatogramExtractor.extractFragmentChromatograms(parameters.getFragmentTolerance(), entry.getLocalizationIons(), downcastedSpectra, null, GraphType.boldline, true, false);
					fragmentTraceMap.putAll(targetMap);
					ArrayList<XYTrace> traces=new ArrayList<XYTrace>(fragmentTraceMap.values());
					
					ChartPanel fragmentChart=Charter.getChart("Retention Time (min)", "Intensity", true, traces.toArray(new XYTrace[traces.size()]));
					
					String annotation = entry.getPeptideAnnotation().getNonDirectionalPeptideAnnotation();
					tabs.add(annotation+" ("+(Math.round(entry.getLocalizationScore()*10.0f)/10.0f)+")", fragmentChart);
					Logger.logLine("Finished reading peptide "+entry.getSpectrumName()+", "+entry.getNumberOfModifiableResidues()+" residues, "+entry.getLocalizationIons().length+" ions");
				}
				chartSplit.setBottomComponent(tabs);

			} catch (Exception e) {
				JOptionPane.showMessageDialog(LocalizationResultsBrowserPanel.this, "Sorry, there was a problem reading the precursor window that contains ["+entries.get(0).getPrecursorMZ()+"]: "+e.getMessage(), "Error Reading DIA File",
						JOptionPane.ERROR_MESSAGE);
				e.printStackTrace();
			}
		}
		tableDataSplit.setDividerLocation(location);
		chartSplit.setDividerLocation(locationRaw);
	}

	private class StatusColumnCellRenderer extends DefaultTableCellRenderer {
		private static final long serialVersionUID = 1L;

		@Override
		public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {
			JLabel l=(JLabel)super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);
			
	        LocalizedLibraryEntryTableModel tableModel = (LocalizedLibraryEntryTableModel) table.getModel();
			if (tableModel.flagRow(row)) {
				l.setBackground(Color.pink);
			} else {
				l.setBackground(Color.white);
			}
			l.setForeground(Color.black);

			return l;

		}
	}
}
