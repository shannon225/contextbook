package edu.washington.gs.maccoss.encyclopedia.gui.framework.pecan;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Optional;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JComboBox;
import javax.swing.JEditorPane;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;

import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorVersion;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.ParametersPanelInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.library.EncyclopediaParametersPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessorTableModel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LabeledComponent;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingWorkerJob;
import edu.washington.gs.maccoss.encyclopedia.jobs.PecanJob;
import edu.washington.gs.maccoss.encyclopedia.jobs.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

public class PecanParametersPanel extends JPanel implements ParametersPanelInterface {
	private static final long serialVersionUID=1L;
	private static final int numberOfCores=Runtime.getRuntime().availableProcessors();
	public static final ImageIcon smallimage=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/pecan_small_icon.png"));
	public static final ImageIcon image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/pecan_icon.png"));
	private static final String programShortDescription="Walnut: PeCAn-based Peptide Detection Directly from Data-Independent Acquisition (DIA) MS/MS Data";
	private static final String programName="Walnut";

	public static final String copy="<html><b><p style=\"font-size:16px; font-family: Helvetica, sans-serif\">Walnut: PeCAn-based Peptide Detection Directly from Data-Independent Acquisition (DIA) MS/MS Data<br></p></b>"
			+ "<p style=\"font-size:10px; font-family: Helvetica, sans-serif\">Walnut uses PeCAn-style scoring to extract peptide fragmentation chromatograms from MZML files, assign peaks, and calculate various peak features. These features are interpreted by Percolator to identify peptides.";

	private static final String[] NUMBER_OF_EXTRA_DECOY_ITEMS=new String[] {"Normal Target/Decoy", "+10% Extra Decoys", "+20% Extra Decoys", "+50% Extra Decoys", "+100% Extra Decoys (2x Time)"};
	private static final float[] NUMBER_OF_EXTRA_DECOY_VALUES=new float[] {0.0f, 0.1f, 0.2f, 0.5f, 1.0f};
	
	private final FileChooserPanel backgroundFasta;
	private final FileChooserPanel targetFasta;
	private final JComboBox<String> fixed=new JComboBox<String>(new String[] {"C+57 (Carbamidomethyl)", "C+58 (Carboxymethyl)", "C+46 (MMTS)", "C+125 (NEM)", "None"});
	private final JComboBox<String> fragType=new JComboBox<String>(new String[] {FragmentationType.toName(FragmentationType.CID), FragmentationType.toName(FragmentationType.HCD), //FragmentationType.toName(FragmentationType.SILAC), 
			FragmentationType.toName(FragmentationType.ETD)});
	private final JComboBox<PercolatorVersion> percolatorVersion=new JComboBox<PercolatorVersion>(PercolatorVersion.VALID_VERSIONS);
	
	private final JFormattedTextField precursorWindowWidth=new JFormattedTextField(NumberFormat.getNumberInstance());

	private final SpinnerModel minCharge=new SpinnerNumberModel(2, 1, 2, 1);
	private final SpinnerModel maxCharge=new SpinnerNumberModel(3, 2, 4, 1);
	private final SpinnerModel maxMissedCleavage=new SpinnerNumberModel(1, 0, 3, 1);
	private final SpinnerModel numberOfQuantitativeIons=new SpinnerNumberModel(5, 1, 100, 1);
	private final SpinnerModel minNumOfQuantitativeIons=new SpinnerNumberModel(3, 0, 100, 1);
	private final JComboBox<String> numberOfExtraDecoyLibraries=new JComboBox<String>(NUMBER_OF_EXTRA_DECOY_ITEMS);
	private final JTextField additionalCommandLineOptions=new JTextField();

	private final SearchPanel searchPanel;
	public PecanParametersPanel(SearchPanel searchPanel) {
		super(new BorderLayout());
		this.searchPanel=searchPanel;

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
		options.add(new LabeledComponent("Target/Decoy Approach", numberOfExtraDecoyLibraries));
		options.add(new LabeledComponent("Fixed", fixed));
		options.add(new LabeledComponent("Maximum Missed Cleavage", new JSpinner(maxMissedCleavage)));
		options.add(new LabeledComponent("Number of Quantitative Ions", new JSpinner(numberOfQuantitativeIons)));
		options.add(new LabeledComponent("Minimum Number of Quantitative Ions", new JSpinner(minNumOfQuantitativeIons)));

		JPanel chargeRange=new JPanel(new FlowLayout());
		chargeRange.setOpaque(true);
		chargeRange.setBackground(Color.white);
		chargeRange.add(new JSpinner(minCharge));
		chargeRange.add(new JLabel("<html><p style=\"font-size:10px; font-family: Helvetica, sans-serif\"> to "));
		chargeRange.add(new JSpinner(maxCharge));
		options.add(new LabeledComponent("Charge range", chargeRange));

		options.add(new LabeledComponent("Additonal Command Line Options", additionalCommandLineOptions));

		this.add(options, BorderLayout.CENTER);
	}

	public ProgramType getProgram() {
		return ProgramType.PecanPie;
	}
	
	public String getProgramShortDescription() {
		return programShortDescription;
	}
	
	public ImageIcon getSmallImage() {
		return smallimage;
	}
	
	public ImageIcon getImage() {
		return image;
	}
	
	@Override
	public void askForSetupFile() {
		backgroundFasta.askForFiles();
	}

	@Override
	public File getBackgroundFastaFile() {
		return backgroundFasta.getFile();
	}
	
	@Override
	public boolean canCreateChromatogramLibraries() {
		return true;
	}
	
	public Optional<String> canLoadData() {
		if (backgroundFasta.getFile()==null) {
			return Optional.of("Please load a background FASTA file first!");
		} else if (targetFasta.getFile()==null) {
			return Optional.of("Please load a target FASTA file first!");		
		}
		return Optional.empty();
	}
	
	public void getJob(File diaFile, JobProcessorTableModel model) {
		PecanSearchParameters parameters=getParameters();
		File fastaFile=backgroundFasta.getFile();
		if (fastaFile==null) return;
		File targetFile=targetFasta.getFile();
		if (targetFile==null) return;
		SearchJob job=getJob(diaFile, fastaFile, targetFile, parameters);

		if (job!=null) {
			model.addJob(job);
		}
	}

	static SearchJob getJob(File diaFile, File fastaFile, File targetFile, PecanSearchParameters parameters) {
		
		ArrayList<FastaPeptideEntry> targets=null;
		if (targetFile!=null&&!targetFile.equals(fastaFile)) {
			Logger.logLine("Reading targets from ["+targetFile.getName()+"]");
			targets=new ArrayList<FastaPeptideEntry>();
			
			ArrayList<FastaEntryInterface> targetProteins=FastaReader.readFasta(targetFile, parameters);
			for (FastaEntryInterface entry : targetProteins) {
				ArrayList<FastaPeptideEntry> peptides=parameters.getEnzyme().digestProtein(entry, parameters.getMinPeptideLength(), parameters.getMaxPeptideLength(), parameters.getMaxMissedCleavages(), parameters.getAAConstants(), false);
				targets.addAll(peptides);
			}
		}
		
		PercolatorExecutionData percolatorFiles=PecanJobData.getPercolatorExecutionData(diaFile, fastaFile, parameters);
		PecanScoringFactory factory=new PecanOneScoringFactory(parameters, percolatorFiles.getInputTSV());
		return new PecanJob(new PecanJobData(Optional.ofNullable(targets), diaFile, fastaFile, percolatorFiles, factory));
	}

	public PecanSearchParameters getParameters() {
		DataAcquisitionType dataAcquisitionType=DataAcquisitionType.DIA;
		DigestionEnzyme digestionEnzyme=searchPanel.getEnzyme();
		AminoAcidConstants aaConstants=AminoAcidConstants.getConstants((String)fixed.getSelectedItem(), new ModificationMassMap());
		FragmentationType fragmentation=FragmentationType.getFragmentationType((String)fragType.getSelectedItem());
		MassTolerance precursorPPMValue=searchPanel.getInstrument().getPrecursorTolerance();
		MassTolerance fragmentPPMValue=searchPanel.getInstrument().getFragmentTolerance();
		byte minChargeValue=((Number)minCharge.getValue()).byteValue();
		byte maxChargeValue=((Number)maxCharge.getValue()).byteValue();
		byte maxMissedCleavageValue=((Number)maxMissedCleavage.getValue()).byteValue();
		Number value=(Number)precursorWindowWidth.getValue();
		float precursorWindowWidthValue=value==null?-1.0f:value.floatValue();
		int numberOfJobsValue=searchPanel.getNumberOfJobs();
		int numberOfQuantitativeIonsValue=((Integer)numberOfQuantitativeIons.getValue());
		int minNumOfQuantitativeIonsValue=((Integer)minNumOfQuantitativeIons.getValue());
		float numberOfExtraDecoyLibrariesValue=NUMBER_OF_EXTRA_DECOY_VALUES[((Integer)numberOfExtraDecoyLibraries.getSelectedIndex())];
		float percolatorThresholdValue=0.01f;

		PercolatorVersion percolator=(PercolatorVersion)percolatorVersion.getSelectedItem();
		PecanSearchParameters parameters=new PecanSearchParameters(
				aaConstants,
				fragmentation,
				precursorPPMValue,
				fragmentPPMValue,
				digestionEnzyme,
				percolator,
				percolatorThresholdValue,
				percolatorThresholdValue,
				PercolatorExecutor.DEFAULT_TRAINING_SET_SIZE,
				PercolatorExecutor.DEFAULT_TRAINING_THRESHOLD,
				PercolatorExecutor.DEFAULT_TRAINING_ITERATIONS,
				maxMissedCleavageValue,
				minChargeValue,
				maxChargeValue,
				dataAcquisitionType,
				precursorWindowWidthValue,
				numberOfJobsValue,
				numberOfQuantitativeIonsValue,
				minNumOfQuantitativeIonsValue,
				-1,
				0.0f,
				numberOfExtraDecoyLibrariesValue,
				true,
				true,
				false,
				searchPanel.getInstrument()
		);

		String cmds=additionalCommandLineOptions.getText();
		HashMap<String, String> params=parameters.toParameterMap();
		params.putAll(CommandLineParser.parseArguments(cmds.split(" ")));
		parameters=PecanParameterParser.parseParameters(params);
		
		return parameters;
	}
	
	public void setParameters(PecanSearchParameters params, String fastaFileName, String targetFileName) {
		if (fastaFileName!=null) {
			File fastaFile=new File(fastaFileName);
			if (fastaFile.exists()) backgroundFasta.update(fastaFile);
		}
		if (targetFileName!=null) {
			File targetFile=new File(targetFileName);
			if (targetFile.exists()) targetFasta.update(targetFile);
		}
		
		fixed.setSelectedItem(AminoAcidConstants.toName(params.getAAConstants()));
		fragType.setSelectedItem(FragmentationType.toName(params.getFragType()));
		
		minCharge.setValue(new Integer(params.getMinCharge()));
		maxCharge.setValue(new Integer(params.getMaxCharge()));
		maxMissedCleavage.setValue(params.getMaxMissedCleavages());
		if (params.getPrecursorWindowSize()>0) {
			precursorWindowWidth.setValue(params.getPrecursorWindowSize());
		} else {
			precursorWindowWidth.setValue(-1);
		}
		int index=Arrays.binarySearch(NUMBER_OF_EXTRA_DECOY_VALUES, params.getNumberOfExtraDecoyLibrariesSearched());
		if (index>=0) {
			numberOfExtraDecoyLibraries.setSelectedIndex(index);
		}
		numberOfQuantitativeIons.setValue(params.getNumberOfQuantitativePeaks());
		minNumOfQuantitativeIons.setValue(params.getMinNumOfQuantitativePeaks());
		percolatorVersion.setSelectedItem(params.getPercolatorVersionNumber());
	}
	
	@Override
	public void savePreferences() {
		try {
			getParameters().savePreferences(backgroundFasta.getFile(), targetFasta.getFile());
		} catch (Exception e) {
			Logger.errorLine("Error writing parameters to disk!");
			Logger.errorException(e);
		}
	}
}
