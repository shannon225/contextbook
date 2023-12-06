package edu.washington.gs.maccoss.encyclopedia.gui.dia.interactive;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Optional;
import java.util.zip.DataFormatException;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
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
import org.jfree.chart.panel.CrosshairOverlay;
import org.jfree.chart.plot.Crosshair;
import org.jfree.chart.plot.XYPlot;
import org.jfree.ui.RectangleEdge;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.dia.FragmentIonConsistencyCharter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingWorkerProgress;
import edu.washington.gs.maccoss.encyclopedia.gui.massspec.ChromatogramCharter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.EditableXYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.ChromatogramExtractor;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.BackgroundSubtractionFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SkylineSGFilter;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

/**
 * click right to approve chromatogram, click left to flag as bad
 * @author searleb
 *
 */
public class ChromatogrindrPanel extends JPanel {
	private static final long serialVersionUID=1L;
	
	private final FileChooserPanel diaFileChooser;
	private final FileChooserPanel libraryFileChooser;
	private final JSplitPane mainSplit=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	private final JSplitPane tableSplit=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	
	private final JTable peptideTable;
	private final TableRowSorter<TableModel> rowSorter;
	private final PeptidePrecursorTableModel peptideModel;
	private final JTextField jtfFilter;
	private final JCheckBox jtfNotFilter=new JCheckBox("NOT");

	private LibraryInterface reference=null;
	private StripeFileInterface dia=null;
	private final SearchParameters parameters;

	
	public static void main(String[] args) {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		final JFrame dialog=new JFrame("Multi ELIB/DIA Detection Browser");

		JMenuBar bar=new JMenuBar();
		JMenu fileMenu=new JMenu("File");
		fileMenu.setMnemonic(KeyEvent.VK_F);
		bar.add(fileMenu);
		
		final ChromatogrindrPanel browser=new ChromatogrindrPanel(params);
		JMenuItem openElib=new JMenuItem("Open Raw Data File...");
		openElib.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				browser.askForData();
			}
		});
		openElib.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		fileMenu.add(openElib);
		
		dialog.setJMenuBar(bar);
		
		dialog.getContentPane().add(browser, BorderLayout.CENTER);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		dialog.pack(); 
		dialog.setSize(1900, 1030);
		dialog.setVisible(true);
	}

	public ChromatogrindrPanel(SearchParameters parameters) {
		super(new BorderLayout());
		
		this.parameters=parameters;
		
		peptideModel=new PeptidePrecursorTableModel();
		peptideTable=new JTable(peptideModel) {
			private static final long serialVersionUID = 1L;

			@Override
			public Object getValueAt(int row, int column) {
				if (column==0) return row+1;
				return super.getValueAt(row, column);
			}
		};
		rowSorter=new TableRowSorter<TableModel>(peptideTable.getModel());
		peptideTable.setRowSorter(rowSorter);

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

		libraryFileChooser=new FileChooserPanel(null, "Reference", new SimpleFilenameFilter(LibraryFile.DLIB, LibraryFile.ELIB), false) {
			private static final long serialVersionUID=1L;

			@Override
			public void update(File... filenames) {
				super.update(filenames);
				if (filenames!=null&&filenames.length>0&&filenames[0]!=null) {
					updateLibrary(filenames[0]);
				}
			}
		};

		diaFileChooser=new FileChooserPanel(null, "Raw file", StripeFileGenerator.getFilenameFilter(), true) {
			private static final long serialVersionUID=1L;

			@Override
			public void update(File... filenames) {
				super.update(filenames);
				if (filenames!=null&&filenames.length>0&&filenames[0]!=null) {
					updateRaw(filenames[0]);
				}
			}
		};
		

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(diaFileChooser);
		options.add(libraryFileChooser);
		
		JPanel buttons=new JPanel(new FlowLayout());
		options.add(buttons);
		
		JButton copyButton = new JButton("Copy");
		copyButton.addActionListener(new ActionListener() {	
			@Override
			public void actionPerformed(ActionEvent e) {
				String copyString=peptideModel.copy();
				StringSelection stringSelection = new StringSelection(copyString);
				Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
				clipboard.setContents(stringSelection, null);
			}
		});
		buttons.add(copyButton);
		
		JButton pasteButton=new JButton("Paste");
		pasteButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
				try {
					String clip=(String)clipboard.getData(DataFlavor.stringFlavor);
					peptideModel.paste(clip);
					
				} catch (IOException | UnsupportedFlavorException ex) {
					Logger.errorLine("Error reading clipboard!");
					Logger.errorException(ex);
				}
			}
		});
		buttons.add(pasteButton);
		
		JPanel searchPanel=new JPanel(new BorderLayout());
		searchPanel.add(new JLabel("Search:"), BorderLayout.WEST);
		searchPanel.add(jtfFilter, BorderLayout.CENTER);
		searchPanel.add(jtfNotFilter, BorderLayout.EAST);

		JPanel left=new JPanel(new BorderLayout());
		left.add(options, BorderLayout.NORTH);
		left.add(new JScrollPane(peptideTable, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED), BorderLayout.CENTER);
		left.add(searchPanel, BorderLayout.SOUTH);
		
		mainSplit.setLeftComponent(left);
		mainSplit.setRightComponent(new JLabel("Select a peptide!"));
		
		setLayout(new BorderLayout());
		add(mainSplit, BorderLayout.CENTER);
	}

	public void askForData() {
		diaFileChooser.askForFiles();
	}

	private void updateFilter() {
		String text=jtfFilter.getText();

		if (text.trim().length()==0) {
			rowSorter.setRowFilter(null);
		} else if (jtfNotFilter.isSelected()) {
			rowSorter.setRowFilter(RowFilter.notFilter(RowFilter.regexFilter("(?i)"+text)));
		} else {
			rowSorter.setRowFilter(RowFilter.regexFilter("(?i)"+text));
		}
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
				reference=t;
				updateToSelectedPeptide();
			}
		};
		worker.execute();
	}

	public void updateRaw(final File f) {
		SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame)SwingUtilities.getWindowAncestor(this), "Please wait...", "Reading Raw File") {
			@Override
			protected Nothing doInBackgroundForReal() throws Exception {
				
				dia=StripeFileGenerator.getFile(f, parameters, true);

				Logger.logLine("Read "+dia.getOriginalFileName()+", ("+dia.getRanges().size()+" total windows)");
				return Nothing.NOTHING;
			}
			@Override
			protected void doneForReal(Nothing t) {
				updateToSelectedPeptide();
			}
		};
		worker.execute();
	}
	
	public void updateToSelectedPeptide() {
		int[] selection=peptideTable.getSelectedRows();
		if (selection.length<=0) return;
		
		InteractivePeptidePrecursor entry=peptideModel.getSelectedRow(peptideTable.convertRowIndexToModel(selection[0]));
		resetPeptide(entry);
	}
	public void resetPeptide(final InteractivePeptidePrecursor entry) {
		int location=mainSplit.getDividerLocation();
		if (location<=5) {
			location=200;
		}
		
		JPanel right=new JPanel(new GridLayout(0, 2));
		right.setBackground(Color.WHITE);
		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		FragmentIon[] primaryIonObjects=model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), false);
		Logger.logLine("Graphing "+entry.getPeptideModSeq()+" ("+primaryIonObjects.length+")"+"...");
		
		try {
			
			ArrayList<XYTrace> fragmentTraces=new ArrayList<>();
			ArrayList<XYTrace> precursorTraces=new ArrayList<>();
			
			float rtInSec=entry.getRetentionTimeInSec();
			float minRTInSec = rtInSec-RT_EXTRACTION_MARGIN_IN_SEC;
			float maxRTInSec = rtInSec+RT_EXTRACTION_MARGIN_IN_SEC;

			Range rtRange=entry.getRTRange();
			if (rtRange!=null) {
				minRTInSec=rtRange.getStart();
				maxRTInSec=rtRange.getStop();
			}
			
			// get precursor traces
			ArrayList<PrecursorScan> precursors=dia.getPrecursors(minRTInSec, maxRTInSec);
			ArrayList<PrecursorScan> trimmedPrecursors=new ArrayList<>();
			for (PrecursorScan spectrum : precursors) {
				if (entry.getPrecursorMZ()>spectrum.getIsolationWindowLower()&&entry.getPrecursorMZ()<spectrum.getIsolationWindowUpper()) {
					trimmedPrecursors.add(spectrum);
				}
			}
			precursors=trimmedPrecursors;
			XYTraceInterface[] traceArray=ChromatogramExtractor.extractPrecursorChromatograms(parameters.getPrecursorTolerance(), entry.getPrecursorMZ(), entry.getPrecursorCharge(), precursors, true, false);
			for (int i = 0; i < traceArray.length; i++) {
				if (traceArray[i] instanceof XYTrace) {
					precursorTraces.add((XYTrace)traceArray[i]);
				}
			}
			
			// get fragment traces
			ArrayList<FragmentScan> scans=dia.getStripes(entry.getPrecursorMZ(), minRTInSec, maxRTInSec, false);
			double[][] allMasses=new double[scans.size()][];
			float[][] allIntensities=new float[scans.size()][];
			float[] retentionTimes=new float[scans.size()];
			for (int i=0; i<scans.size(); i++) {
				FragmentScan scan=scans.get(i);
				Pair<double[], float[]> results=extract(scan, primaryIonObjects);
				double[] masses=results.x;
				float[] intensities=results.y;
				
				allMasses[i]=masses;
				allIntensities[i]=intensities;
				retentionTimes[i]=scan.getScanStartTime();
			}

			int movingAverageLength=8; // expected points across the peak
			float[][] chromatograms=General.transposeMatrix(allIntensities);
			ArrayList<float[]> chromatogramList=new ArrayList<float[]>();
			ArrayList<FragmentIon> foundIons=new ArrayList<>();
			for (int j = 0; j < chromatograms.length; j++) {
				if (General.sum(chromatograms[j])>0.0f) {
					chromatograms[j]=SkylineSGFilter.paddedSavitzkyGolaySmooth(chromatograms[j]);
					if (parameters.isSubtractBackground()) {
						chromatograms[j]=BackgroundSubtractionFilter.backgroundSubtractMovingMedian(chromatograms[j], movingAverageLength*10);
					}
					chromatogramList.add(chromatograms[j]);
					foundIons.add(primaryIonObjects[j]);
				}
			}
			primaryIonObjects=foundIons.toArray(new FragmentIon[0]);
			
			TransitionRefinementData data=TransitionRefiner.identifyTransitions(entry.getPeptideModSeq(), entry.getPrecursorCharge(), entry.getRetentionTimeInSec(), 
					primaryIonObjects, chromatogramList, retentionTimes, false, parameters);
			fragmentTraces=getTraces(data.getChromatograms(), data.getCorrelationArray(), retentionTimes, data.getRange());
			
			double globalMaxYFragment=0.0;
			double globalMaxYPrecursor=0.0;
			for (XYTrace xyTrace : precursorTraces) {
				if (xyTrace.getType()==GraphType.line) {
					globalMaxYPrecursor=Math.max(globalMaxYPrecursor, xyTrace.getMaxY());
				}
			}
			for (XYTrace xyTrace : fragmentTraces) {
				if (xyTrace.getType()==GraphType.boldline) {
					globalMaxYFragment=Math.max(globalMaxYFragment, xyTrace.getMaxY());
				}
			}
			
			final ChartPanel chartPanel=ChromatogramCharter.createChart(Optional.ofNullable(precursorTraces), Optional.ofNullable(fragmentTraces), globalMaxYPrecursor, globalMaxYFragment);
			chartPanel.setMouseZoomable(false, false);

			CrosshairOverlay crosshairOverlay = new CrosshairOverlay();
			final Crosshair xCrosshair=new Crosshair(Double.NaN, Color.GRAY, new BasicStroke(0f));
			final EditableXYPoint zoomPoint=new EditableXYPoint();
	        xCrosshair.setLabelVisible(true);
	        crosshairOverlay.addDomainCrosshair(xCrosshair);
	        chartPanel.addOverlay(crosshairOverlay);
	        chartPanel.mouseDragged(null);
	        
	        chartPanel.addMouseMotionListener(new MouseMotionListener() {
				@Override
				public void mouseMoved(MouseEvent e) {
		            Rectangle2D area=chartPanel.getScreenDataArea(e.getX(), e.getY());
			        if (area!=null) {
			        XYPlot plot = (XYPlot) chartPanel.getChart().getPlot();
			        ValueAxis xAxis = plot.getDomainAxis();
			        	double x = xAxis.java2DToValue(e.getX(), area, RectangleEdge.BOTTOM);
			        	xCrosshair.setValue(x);
			        }
				}
				
				@Override
				public void mouseDragged(MouseEvent e) {
		            Rectangle2D area=chartPanel.getScreenDataArea(e.getX(), e.getY());
			        if (area!=null) {
				        XYPlot plot = (XYPlot) chartPanel.getChart().getPlot();
				        ValueAxis xAxis = plot.getDomainAxis();
				        double x = xAxis.java2DToValue(e.getX(), area, RectangleEdge.BOTTOM);
				        xCrosshair.setValue(x);
				        
				        double prevX=zoomPoint.getX();
				        double prevY=zoomPoint.getY();
				        if (!Double.isNaN(prevX)&&!Double.isNaN(prevY)) {
					        Line2D zoomLine=new Line2D.Double(prevX, prevY, e.getX(), prevY);
	
					        Graphics2D g2 = (Graphics2D) chartPanel.getGraphics();
					        g2.setPaint(Color.gray);
					        g2.draw(zoomLine);
				        }
			        }
				}
			});
	        
	        chartPanel.addMouseListener(new MouseListener() {
				@Override
				public void mouseReleased(MouseEvent e) {
					double prevX=zoomPoint.getX();
					if (!Double.isNaN(prevX)) {
			            Rectangle2D area=chartPanel.getScreenDataArea(e.getX(), e.getY());
		            	double x=Math.max(area.getMinX(), Math.min(e.getX(), area.getMaxX()));
		            	
		            	double first=Math.min(x, prevX);
		            	double second=Math.max(x, prevX);
		            	
				        XYPlot plot = (XYPlot) chartPanel.getChart().getPlot();
				        ValueAxis xAxis = plot.getDomainAxis();
				        double plotX1 = xAxis.java2DToValue(first, area, RectangleEdge.BOTTOM);
				        double plotX2 = xAxis.java2DToValue(second, area, RectangleEdge.BOTTOM);
				        
				        entry.setRtRangeInSecs(new Range(plotX1*60.0, plotX2*60.0));
				        resetPeptide(entry);
			        }
				}
				
				@Override
				public void mousePressed(MouseEvent e) {
		            Rectangle2D area=chartPanel.getScreenDataArea(e.getX(), e.getY());
		            if (area!=null) {
		            	zoomPoint.setX(Math.max(area.getMinX(), Math.min(e.getX(), area.getMaxX())));
		            	zoomPoint.setY(Math.max(area.getMinY(), Math.min(e.getY(), area.getMaxY())));
		            }
		            else {
		            	zoomPoint.setX(null);
		            	zoomPoint.setY(null);
		            }
				}
				
				@Override
				public void mouseExited(MouseEvent e) {
				}
				
				@Override
				public void mouseEntered(MouseEvent e) {
				}
				
				@Override
				public void mouseClicked(MouseEvent e) {
				}
			});
	        
			right.add(chartPanel);
			
			if (reference!=null) {
				ArrayList<LibraryEntry> references=reference.getEntries(entry.getPeptideModSeq(), entry.getPrecursorCharge(), false);
				if (references.size()==0) {
					references=reference.getEntries(entry.getLegacyPeptideModSeq(), entry.getPrecursorCharge(), false);
				}
				
				if (references.size()>0) {
					LibraryEntry ref=references.get(0);
					LibraryEntry acq=data.getEntry(ref, parameters);
	
					LibraryEntry butterfly=FragmentIonConsistencyCharter.getButterfly(acq, ref);
					ChartPanel chartPanelButterfly = Charter.getChart(new AnnotatedLibraryEntry(butterfly, parameters, true));

					right.add(chartPanelButterfly);
				}
			}
	
			mainSplit.setRightComponent(right);

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
		
		mainSplit.setDividerLocation(location);
	}
	
	private static final float RT_EXTRACTION_MARGIN_IN_SEC=45f;

	public Pair<double[], float[]> extract(Spectrum spectrum, FragmentIon[] ions) {
		MassTolerance acquiredTolerance=parameters.getFragmentTolerance();
		
		double[] acquiredMasses=spectrum.getMassArray();
		float[] acquiredIntensities=spectrum.getIntensityArray();

		TDoubleArrayList actualTargetMasses=new TDoubleArrayList();
		TFloatArrayList actualTargetIntensities=new TFloatArrayList();
		for (int i = 0; i < ions.length; i++) {
			FragmentIon target=ions[i];
		
			int[] indicies=acquiredTolerance.getIndicies(acquiredMasses, target.getMass());
			float intensity=0.0f;
			float bestPeakIntensity=0.0f;
			double bestPeakMass=0.0;
			
			for (int j=0; j<indicies.length; j++) {
				intensity+=acquiredIntensities[indicies[j]];
				
				if (acquiredIntensities[indicies[j]]>bestPeakIntensity) {
					bestPeakIntensity=acquiredIntensities[indicies[j]];
					bestPeakMass=acquiredMasses[indicies[j]];
				}
			}
			actualTargetIntensities.add(intensity);
			actualTargetMasses.add(bestPeakMass);
		}

		double[] actualTargetMassesRaw=actualTargetMasses.toArray();
		float[] actualTargetIntensitiesRaw=actualTargetIntensities.toArray();
		
		return new Pair<double[], float[]>(actualTargetMassesRaw, actualTargetIntensitiesRaw);
	}

	private static ArrayList<XYTrace> getTraces(ArrayList<float[]> chromatograms, float[] correlationArray, float[] rts, Range rtRange) {
		ArrayList<XYTrace> xytraces=new ArrayList<XYTrace>();
		for (int i=0; i<chromatograms.size(); i++) {
			float[] fs=chromatograms.get(i);
			
			Color c;
			GraphType graphtype;
			float thickness;
			if (correlationArray[i]>TransitionRefiner.quantitativeCorrelationThreshold) {
				c=new Color(0, 205, 0);
				graphtype=GraphType.boldline;
				thickness=3.0f;
			} else if (correlationArray[i]>TransitionRefiner.identificationCorrelationThreshold) {
				c=new Color(255, 215, 0);
				graphtype=GraphType.boldline;
				thickness=3.0f;
			} else if (correlationArray[i]==0.0f) {
				c=Color.gray;
				graphtype=GraphType.line;
				thickness=1.0f;
			} else {
				c=Color.red;
				graphtype=GraphType.dashedline;
				thickness=3.0f;
			}

			xytraces.add(TransitionRefiner.toXYTrace(fs, rts, ""+i, c, rtRange, graphtype, thickness));

			if (rtRange!=null) {
				xytraces.add(TransitionRefiner.toXYTrace(fs, rts, ""+i, Color.gray, null, graphtype, thickness));
			}
		}
		return xytraces;
	}

}
