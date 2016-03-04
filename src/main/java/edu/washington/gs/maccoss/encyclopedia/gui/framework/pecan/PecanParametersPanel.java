package edu.washington.gs.maccoss.encyclopedia.gui.framework.pecan;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JComboBox;
import javax.swing.JEditorPane;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessorTableModel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LabeledComponent;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingJob;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;

public class PecanParametersPanel extends JPanel {
	private static final long serialVersionUID=1L;
	private static final int numberOfCores=Runtime.getRuntime().availableProcessors();
	public static final ImageIcon image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/pecan.png"));
	public static final String copy="<html><b><p style=\"font-size:16px; font-family: Helvetica, sans-serif\">PECAN: Peptide Detection Directly from Data-Independent Acquisition (DIA) MS/MS Data<br></p></b>"
			+ "<p style=\"font-size:10px; font-family: Helvetica, sans-serif\">PECAN extracts peptide fragmentation chromatograms from MZML files, assigns peaks, and calculates various peak features. These features are interpreted by Percolator to identify peptides.";
	
	
	private final FileChooserPanel backgroundFasta;
	private final FileChooserPanel targetFasta;
	private final JComboBox<String> overlap=new JComboBox<String>(new String[] {"Overlapped", "Not Overlapped"});
	private final JComboBox<String> enzyme=new JComboBox<String>(new String[] {"Trypsin", "Lys-C", "Lys-N", "Arg-C", "CNBr", "Chymotrypsin", "PepsinA"});
	private final JComboBox<String> fixed=new JComboBox<String>(new String[] {"C+57 (Carbamidomethyl)", "C+58 (Carboxymethyl)", "C+46 (MMTS)", "None"});
	private final JComboBox<String> fragType=new JComboBox<String>(new String[] {"HCD (Y-Only)", "CID (B/Y)", "ETD (C/Z/Z+1)"});

	private final SpinnerModel precursorPPM = new SpinnerNumberModel(10, 1, 1000, 1);
	private final SpinnerModel fragmentPPM = new SpinnerNumberModel(10, 1, 1000, 1);
	private final SpinnerModel minCharge = new SpinnerNumberModel(2, 1, 2, 1);
	private final SpinnerModel maxCharge = new SpinnerNumberModel(3, 2, 4, 1);
	private final SpinnerModel maxMissedCleavage = new SpinnerNumberModel(1, 0, 3, 1);
	private final SpinnerModel numberOfJobs = new SpinnerNumberModel(numberOfCores, 1, numberOfCores, 1);

	public PecanParametersPanel() {
		super(new BorderLayout());

		JPanel top=new JPanel(new BorderLayout());
		top.add(new JLabel(image), BorderLayout.WEST);
		JEditorPane editor=new JEditorPane("text/html", copy);
		editor.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
		top.add(editor, BorderLayout.CENTER);
		top.setBackground(Color.white);
		this.add(top, BorderLayout.NORTH);
		
		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(new LabeledComponent("<p style=\"font-size:12px; font-family: Helvetica, sans-serif\"><b>Parameters", new JLabel()));
		
		backgroundFasta=new FileChooserPanel(null, "Background", new SimpleFilenameFilter(".fas", ".fasta"), true) {
			private static final long serialVersionUID=1L;

			@Override
			public void update(File... filename) {
				super.update(filename);
				if (filename!=null&&filename.length>0&&filename[0]!=null) {
					if (targetFasta.getFile()==null) {
						targetFasta.update(filename);
					}
				}
			}
		};
		options.add(backgroundFasta);
		targetFasta=new FileChooserPanel(null, "Target", new SimpleFilenameFilter(".fas", ".fasta"), true);
		options.add(targetFasta);
		options.add(new LabeledComponent("Precursor Isolation Windows", overlap));
		options.add(new LabeledComponent("Enzyme", enzyme));
		options.add(new LabeledComponent("Fixed", fixed));
		options.add(new LabeledComponent("Fragmentation", fragType));
		options.add(new LabeledComponent("Precursor (PPM)", new JSpinner(precursorPPM)));
		options.add(new LabeledComponent("Fragment (PPM)", new JSpinner(fragmentPPM)));
		options.add(new LabeledComponent("Maximum Missed Cleavage", new JSpinner(maxMissedCleavage)));
		options.add(new LabeledComponent("Number of Cores", new JSpinner(numberOfJobs)));

		JPanel chargeRange=new JPanel(new FlowLayout());
		chargeRange.setOpaque(true);
		chargeRange.setBackground(Color.white);
		chargeRange.add(new JSpinner(minCharge));
		chargeRange.add(new JLabel("<html><p style=\"font-size:10px; font-family: Helvetica, sans-serif\"> to "));
		chargeRange.add(new JSpinner(maxCharge));
		options.add(new LabeledComponent("Charge range", chargeRange));
		

		this.add(options, BorderLayout.CENTER);
	}
	
	public Optional<String> canLoadData() {
		if (backgroundFasta.getFile()==null) {
			return Optional.of("Please load a background FASTA file first!");
		} else if (targetFasta.getFile()==null) {
			return Optional.of("Please load a target FASTA file first!");		
		}
		return Optional.empty();
	}
	
	public SwingJob getJob(File diaFile, JobProcessorTableModel model) {
		PecanSearchParameters parameters=getParameters();
		File fastaFile=backgroundFasta.getFile();
		if (fastaFile==null) return null;
		File targetFile=targetFasta.getFile();
		if (targetFile==null) return null;
		SearchJob job=getJob(diaFile, fastaFile, targetFile, model, parameters);

		if (job!=null) {
			model.addJob(job);
		}
		return job;
	}

	static SearchJob getJob(File diaFile, File fastaFile, File targetFile, JobProcessor processor, PecanSearchParameters parameters) {
		File outputFile=new File(diaFile.getAbsolutePath()+".pecan.txt");
		File featureFile=new File(outputFile.getAbsolutePath()+".features.txt");
		
		ArrayList<FastaEntry> targets=null;
		if (targetFile!=null&&!targetFile.equals(fastaFile)) {
			Logger.logLine("Reading targets from ["+targetFile.getName()+"]");
			targets=new ArrayList<FastaEntry>();
			
			ArrayList<FastaEntry> targetProteins=FastaReader.readFasta(targetFile);
			for (FastaEntry entry : targetProteins) {
				ArrayList<String> peptides=parameters.getEnzyme().digestProtein(entry.getSequence(), parameters.getMinPeptideLength(), parameters.getMaxPeptideLength(), parameters.getMaxMissedCleavages());
				for (String peptide : peptides) {
					FastaEntry pe=entry.getSubEntry(peptide);
					targets.add(pe);
				}
			}
		}
		
		PecanScoringFactory factory=new PecanOneScoringFactory(parameters, featureFile);
		return new PecanJob(processor, new PecanJobData(Optional.ofNullable(targets), diaFile, fastaFile, featureFile, outputFile, factory));
	}

	private PecanSearchParameters getParameters() {
		boolean deconvoluteOverlappingWindows="Overlapped".equalsIgnoreCase((String)overlap.getSelectedItem());
		DigestionEnzyme digestionEnzyme=DigestionEnzyme.getEnzyme((String)enzyme.getSelectedItem());
		AminoAcidConstants aaConstants=AminoAcidConstants.getConstants((String)fixed.getSelectedItem());
		FragmentationType fragmentation=FragmentationType.getFragmentationType((String)fragType.getSelectedItem());
		float precursorPPMValue=((Integer)precursorPPM.getValue()).floatValue();
		float fragmentPPMValue=((Integer)fragmentPPM.getValue()).floatValue();
		byte minChargeValue=((Integer)minCharge.getValue()).byteValue();
		byte maxChargeValue=((Integer)maxCharge.getValue()).byteValue();
		byte maxMissedCleavageValue=((Integer)maxMissedCleavage.getValue()).byteValue();
		int numberOfJobsValue=((Integer)numberOfJobs.getValue());
		PecanSearchParameters parameters=new PecanSearchParameters(aaConstants, fragmentation, new MassTolerance(precursorPPMValue), new MassTolerance(fragmentPPMValue), digestionEnzyme,
				maxMissedCleavageValue, minChargeValue, maxChargeValue, deconvoluteOverlappingWindows, numberOfJobsValue);
		return parameters;
	}
}
