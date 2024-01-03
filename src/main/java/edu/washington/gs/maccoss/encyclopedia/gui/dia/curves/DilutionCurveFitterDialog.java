package edu.washington.gs.maccoss.encyclopedia.gui.dia.curves;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FilenameFilter;
import java.util.HashSet;
import java.util.StringTokenizer;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTextArea;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;

import edu.washington.gs.maccoss.encyclopedia.algorithms.curve.AbstractDilutionCurveFittingParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.curve.DilutionCurveFitter;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LabeledComponent;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;

public class DilutionCurveFitterDialog extends JDialog {
	private static final long serialVersionUID = 1L;
	private static final ImageIcon processingIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/processing_icon.png"));

	public static void launchDilutionCurveFitterDialog(Component root) {
		JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		
		final JDialog dialog=new JDialog(frame, "DIA to PRM Dialog", true);

		final FileChooserPanel libraryFileChooser=new FileChooserPanel(null, "Global Library", new SimpleFilenameFilter(".dlib", ".elib"), true);
		final FileChooserPanel singleFileChooser=new FileChooserPanel(null, "Alignment Search Result", new SimpleFilenameFilter(".dlib", ".elib"), true);
		final FileChooserPanel titrationCurveFileChooser=new FileChooserPanel(null, "Titration Curve File", new SimpleFilenameFilter(".tsv", ".csv", ".txt"), true);
		final FileChooserPanel sampleOrganizationFileChooser=new FileChooserPanel(null, "Sample Organization File", new SimpleFilenameFilter(".tsv", ".csv", ".txt"), true);
		final FileChooserPanel saveDirFileChooser=new FileChooserPanel(null, "Save Directory", new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				if (dir.exists()&&!dir.isDirectory()) {
					return false;
				}
				return true;
			}
		}, true, false);
		
		final JTextArea eliminateTextArea = new JTextArea(25, 40);
		eliminateTextArea.setFont(new Font("Monospaced", Font.PLAIN, 10));
		JScrollPane eliminateScrollPane = new JScrollPane(eliminateTextArea); 
		eliminateScrollPane.setBorder(BorderFactory.createTitledBorder("Eliminated peptides (poor reproducibility, etc):"));
		eliminateTextArea.setEditable(true);
		
		final JTextArea targetTextArea = new JTextArea(25, 40);
		targetTextArea.setFont(new Font("Monospaced", Font.PLAIN, 10));
		JScrollPane targetScrollPane = new JScrollPane(targetTextArea);
		targetScrollPane.setBorder(BorderFactory.createTitledBorder("Target protein accessions (requires exact matches):"));
		targetTextArea.setEditable(true);

		final SpinnerModel maximumPeptidesPerProteinSpinner=new SpinnerNumberModel(3, 1, 100, 1);
		final SpinnerModel maximumAssayDensitySpinner=new SpinnerNumberModel(10, 1, 100, 1);
		final SpinnerModel windowSizeInMinSpinner=new SpinnerNumberModel(5.0, 0.1, 100.0, 1.0);

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(libraryFileChooser);
		options.add(singleFileChooser);
		options.add(new LabeledComponent("Maximum number of peptides per protein", new JSpinner(maximumPeptidesPerProteinSpinner)));
		options.add(new LabeledComponent("Maximum assay density (scans/second)", new JSpinner(maximumAssayDensitySpinner)));
		options.add(new LabeledComponent("Retention time window size (min)", new JSpinner(windowSizeInMinSpinner)));
		
		options.add(targetScrollPane);
		options.add(eliminateScrollPane);

		options.add(saveDirFileChooser);
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchParameters params=SearchParameterParser.getDefaultParametersObject();
				final File saveDir=saveDirFileChooser.getFile();
				final File libraryFile=libraryFileChooser.getFile();
				final File exampleSearchFile=singleFileChooser.getFile();
				final File titrationCurveFile=titrationCurveFileChooser.getFile();
				final File sampleOrganizationFile=sampleOrganizationFileChooser.getFile();

				final String[] targets = getTokens(targetTextArea.getText());
				final String[] eliminatedPeptides = getTokens(eliminateTextArea.getText());
				
				if (saveDir==null) {
					JOptionPane.showMessageDialog(frame, "You must specify a global ELIB or DLIB library file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, processingIcon);
				} else if (targets.length==0) {
					JOptionPane.showMessageDialog(frame, "You must specify at least one target accession number!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, processingIcon);
				} else if (exampleSearchFile==null) {
					JOptionPane.showMessageDialog(frame, "You must specify a retention time alignment file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, processingIcon);
				} else if (titrationCurveFile==null) {
					JOptionPane.showMessageDialog(frame, "You must specify a titration curve file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, processingIcon);
				} else if (sampleOrganizationFile==null) {
					JOptionPane.showMessageDialog(frame, "You must specify a sample organization file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, processingIcon);
				} else {
					dialog.setVisible(false);
					dialog.dispose();
					
					DilutionCurveParameters fittingParams=new DilutionCurveParameters(((Number)maximumPeptidesPerProteinSpinner.getValue()).intValue(), 
							((Number)windowSizeInMinSpinner.getValue()).floatValue(), ((Number)maximumAssayDensitySpinner.getValue()).intValue(), targets, eliminatedPeptides);

					try {
						DilutionCurveFitter.generateAssayFromCurves(params, saveDir, titrationCurveFile, sampleOrganizationFile, libraryFile, exampleSearchFile, fittingParams);
					} catch (Exception ex) {
						JOptionPane.showMessageDialog(frame, "Unexpected error: "+ex.getMessage(), "Unexpected Error!", JOptionPane.ERROR_MESSAGE, processingIcon);
					}
				}
			}

			private String[] getTokens(final String text) {
				final HashSet<String> targets=new HashSet<>();
				if (text!=null&&text.length()>0) {
					StringTokenizer st=new StringTokenizer(text);
					while (st.hasMoreTokens()) {
						targets.add(st.nextToken());
					}
				}
				return targets.toArray(new String[0]);
			}
		});
		buttons.add(okButton);
		JButton cancelButton=new JButton("Cancel");
		cancelButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				dialog.setVisible(false);
				dialog.dispose();
			}
		});
		buttons.add(cancelButton);
		
		JPanel mainpane=new JPanel(new BorderLayout());
		mainpane.add(options, BorderLayout.CENTER);
		mainpane.add(buttons, BorderLayout.SOUTH);
		mainpane.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10), BorderFactory.createTitledBorder("Parameters:")));
		
		dialog.getContentPane().add(mainpane, BorderLayout.CENTER);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		dialog.pack(); 
		dialog.setSize(500, 700);
		dialog.setVisible(true);
	}

	private static class DilutionCurveParameters implements AbstractDilutionCurveFittingParameters {
		private final int numberOfRTAnchors=0;
		private final int maxNumberPeptidesPerProtein;
		private final int targetTotalNumberOfPeptides=300; // remember to subtract off anchors (total is 160 peptides)
		private final float windowInMin; // in minutes!
		private final float minCVForAnchors=0.05f;
		private final float minCVForBadAnchors=0.75f;
		private final int assayMaxDensity;
		private final boolean requireAlignmentRT=false; // turn off for fitting against PRM
		private final boolean useLineNoise=true; // newer versions should set this to "true"
		private final String[] targets; 
		private final String[] eliminatedPeptides; 
		
		
		public DilutionCurveParameters(int maxNumberPeptidesPerProtein, float windowInMin, int assayMaxDensity,
				String[] targets, String[] eliminatedPeptides) {
			super();
			this.maxNumberPeptidesPerProtein = maxNumberPeptidesPerProtein;
			this.windowInMin = windowInMin;
			this.assayMaxDensity = assayMaxDensity;
			this.targets = targets;
			this.eliminatedPeptides = eliminatedPeptides;
		}
		
		public float getWindowInMin() {
			return windowInMin;
		}
		public float getWindowInMin(float rtInSec) {
			return windowInMin;
		}

		public int getNumberOfRTAnchors() {
			return numberOfRTAnchors;
		}

		public int getMaxNumberPeptidesPerProtein() {
			return maxNumberPeptidesPerProtein;
		}

		public int getTargetTotalNumberOfPeptides() {
			return targetTotalNumberOfPeptides;
		}

		public float getMinCVForAnchors() {
			return minCVForAnchors;
		}

		public float getMinCVForBadAnchors() {
			return minCVForBadAnchors;
		}

		public int getAssayMaxDensity() {
			return assayMaxDensity;
		}

		@Override
		public boolean isTargetedProtein(String accession) {
			for (int i = 0; i < targets.length; i++) {
				if (accession.indexOf(targets[i])>=0) return true;
			}
			return false;
		}

		@Override
		public boolean isEliminatedPeptide(String peptideModSeq) {
			for (int i = 0; i < eliminatedPeptides.length; i++) {
				if (peptideModSeq.indexOf(eliminatedPeptides[i])>=0) return true;
			}
			return false;
		}

		public boolean isRequireAlignmentRT() {
			return requireAlignmentRT;
		}

		public boolean isUseLineNoise() {
			return useLineNoise;
		}
	}
}
