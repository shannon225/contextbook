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
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JSplitPane;
import javax.swing.JTextField;
import javax.swing.JToolBar;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.UIManager;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ExpectedFragmentationScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.FragmentationScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.FragmentationTraceTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneFragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanRawScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScoringTask;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SkylineSGFilter;

public class PeptideExtractingBrowserPanel extends JPanel {
	private static final long serialVersionUID=1L;
	
	private final PecanRawScorer scorer;

	private final SearchParameters parameters;
	private final FileChooserPanel diaFile;
	private final JTextField peptide=new JTextField("YLDGLTAER");
	private final SpinnerModel charge=new SpinnerNumberModel(2, 1, 5, 1);
	private final JSplitPane split=new JSplitPane(JSplitPane.VERTICAL_SPLIT);

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

		JToolBar bar=new JToolBar();
		bar.add(diaFile);
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
		add(split, BorderLayout.CENTER);
	}

	public void resetPeptide(String peptide, int charge) {
		if (peptide==null||peptide.length()==0||dia==null) {
			split.setTopComponent(new JLabel("Select a peptide!"));
			split.setBottomComponent(new JPanel());
		} else {
			Logger.logLine("Parsing peptide...");
			PecanOneFragmentationModel model=new PecanOneFragmentationModel(new FastaPeptideEntry(peptide), parameters.getAAConstants());
			ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
			LibraryEntry entry=model.getUnitSpectrum(dia.getOriginalFileName(), new HashSet<String>(), (byte)charge, 0.0f, parameters);
			entries.add(entry);
			
			try {
				List<Stripe> stripes=dia.getStripes(entry.getPrecursorMZ(), 0.0f, Float.MAX_VALUE, false);
				FragmentationTraceTask task=new FragmentationTraceTask(scorer, FragmentationTraceTask.PLOT_INTENSITIES, entries, stripes, new PrecursorScanMap(new ArrayList<PrecursorScan>()), parameters.getAAConstants());
				HashMap<LibraryEntry, PeptideScoringResult> result=task.call();
				
				ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
				for (Entry<LibraryEntry, PeptideScoringResult> resultEntry : result.entrySet()) {
					FragmentationScoringResult peptideResult=(FragmentationScoringResult)resultEntry.getValue();

					for (XYTrace trace : peptideResult.getFragmentationTraces()) {
						XYTrace sgSmoothed=SkylineSGFilter.paddedSavitzkyGolaySmooth(trace);
						traces.add(sgSmoothed);
					}
				}
				ChartPanel chart=Charter.getChart("Retention Time (min)", "Intensity", true, traces.toArray(new XYTrace[traces.size()]));
				split.setTopComponent(chart);
				

				
				/*BlockingQueue<PeptideScoringResult> ionCountResultsQueue=new LinkedBlockingQueue<PeptideScoringResult>();
				IonCountingScoringTask ionCount=new IonCountingScoringTask(scorer, entries, stripes, 2.5f, new PrecursorScanMap(new ArrayList<PrecursorScan>()), ionCountResultsQueue, parameters);
				ionCount.call();
				
				PeptideScoringResult ionCountResult=ionCountResultsQueue.take();
				XYTraceInterface xytrace=ionCountResult.getTrace();
				Pair<double[], double[]> trace=xytrace.toArrays();
				double[] newx=General.multiply(trace.x, 1.0f/60.0f); // scale to minutes
				XYTraceInterface ionCounttrace=new XYTrace(newx, trace.y, xytrace.getType(), xytrace.getName(), xytrace.getColor(), xytrace.getThickness());
				
				ChartPanel ionCountchart=Charter.getChart("Retention Time (min)", "RawScore", false, ionCounttrace);*/
				BlockingQueue<PeptideScoringResult> resultsQueue=new LinkedBlockingQueue<PeptideScoringResult>();
				XCorDIAOneScoringTask xcorrTask=new XCorDIAOneScoringTask(new XCorDIAOneScorer(parameters, null), entries, stripes, 2.5f, new PrecursorScanMap(new ArrayList<PrecursorScan>()), resultsQueue, parameters);
				xcorrTask.call();
				
				PeptideScoringResult ionCountResult=resultsQueue.take();
				XYTraceInterface xytrace=ionCountResult.getTrace();
				Pair<double[], double[]> trace=xytrace.toArrays();
				double[] newx=General.multiply(trace.x, 1.0f/60.0f); // scale to minutes
				XYTraceInterface ionCounttrace=new XYTrace(newx, trace.y, xytrace.getType(), xytrace.getName(), xytrace.getColor(), xytrace.getThickness());
				
				ChartPanel ionCountchart=Charter.getChart("Retention Time (min)", "XCorr", false, ionCounttrace);

				split.setBottomComponent(ionCountchart);
				
				
			} catch (Exception e) {
				e.printStackTrace();
				JOptionPane.showMessageDialog(PeptideExtractingBrowserPanel.this, "Sorry, there was a problem reading the precursor window that contains ["+entry.getPrecursorMZ()+"]: "+e.getMessage(), "Error Reading DIA File",
						JOptionPane.ERROR_MESSAGE);
			}
			Logger.logLine("Finished reading peptide "+peptide);
		}
	}
}
