package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JSplitPane;
import javax.swing.JTextField;
import javax.swing.JToolBar;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ExpectedFragmentationScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.FragmentationScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.FragmentationTraceTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneFragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanRawScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;

public class DIABrowserPanel extends JPanel {
	private static final long serialVersionUID=1L;
	
	private final PecanRawScorer scorer;

	private final SearchParameters parameters;
	private final FileChooserPanel diaFile;
	private final JTextField peptide=new JTextField("YLDGLTAER");
	private final SpinnerModel charge=new SpinnerNumberModel(2, 1, 5, 1);
	private final JSplitPane split=new JSplitPane(JSplitPane.VERTICAL_SPLIT);

	private StripeFileInterface dia=null;

	public DIABrowserPanel() {
		super(new BorderLayout());
		HashMap<String, String> map=SearchParameterParser.getDefaultParameters();
		map.put("-deconvoluteOverlappingWindows", "true");
		parameters=SearchParameterParser.parseParameters(map);
		scorer=new PecanRawScorer(parameters.getFragmentTolerance(), new ExpectedFragmentationScorer(parameters));

		diaFile=new FileChooserPanel(null, "DIA File", new SimpleFilenameFilter(".dia", ".mzml")) {
			private static final long serialVersionUID=1L;

			@Override
			public void update(File... filename) {
				super.update(filename);
				if (filename!=null&&filename.length>0&&filename[0]!=null) {
					try {
						Logger.logLine("Reading file...");

						dia=MzmlToDIAConverter.getFile(filename[0], parameters);
						Logger.logLine("Finished reading file.");
						resetPeptide(peptide.getText(), (Integer) charge.getValue());
					} catch (Exception e) {
						e.printStackTrace();
						JOptionPane.showMessageDialog(DIABrowserPanel.this, "Sorry, there was a problem reading ["+filename[0].getName()+"]: "+e.getMessage(), "Error Opening DIA File",
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
			PecanOneFragmentationModel model=new PecanOneFragmentationModel(new FastaEntry(peptide), parameters.getAAConstants());
			ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
			LibraryEntry entry=model.getUnitSpectrum((byte)charge, parameters);
			entries.add(entry);
			
			try {
				ArrayList<Stripe> stripes=dia.getStripes(entry.getPrecursorMZ(), 0.0f, Float.MAX_VALUE, false);
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
				split.setBottomComponent(chart);
				
			} catch (Exception e) {
				JOptionPane.showMessageDialog(DIABrowserPanel.this, "Sorry, there was a problem reading the precursor window that contains ["+entry.getPrecursorMZ()+"]: "+e.getMessage(), "Error Reading DIA File",
						JOptionPane.ERROR_MESSAGE);
			}
			Logger.logLine("Finished reading peptide "+peptide);
		}
	}
}
