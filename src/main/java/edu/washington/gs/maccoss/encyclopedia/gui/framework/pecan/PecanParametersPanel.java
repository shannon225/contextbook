package edu.washington.gs.maccoss.encyclopedia.gui.framework.pecan;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
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
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.ParametersPanelInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.library.EncyclopediaParametersPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessorTableModel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LabeledComponent;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingJob;
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
	private final JComboBox<String> acquisition=new JComboBox<String>(new String[] {DataAcquisitionType.toName(DataAcquisitionType.OVERLAPPING_DIA), DataAcquisitionType.toName(DataAcquisitionType.DIA)});
	private final JComboBox<String> enzyme=new JComboBox<String>(new String[] {"Trypsin", "Lys-C", "Lys-N", "Arg-C", "CNBr", "Chymotrypsin", "Pepsin A", "No Enzyme"});
	private final JComboBox<String> fixed=new JComboBox<String>(new String[] {"C+57 (Carbamidomethyl)", "C+58 (Carboxymethyl)", "C+46 (MMTS)", "None"});
	private final JComboBox<String> fragType=new JComboBox<String>(new String[] {FragmentationType.toName(FragmentationType.CID), FragmentationType.toName(FragmentationType.HCD), FragmentationType.toName(FragmentationType.ETD)});
	private final JComboBox<String> percolatorVersion=new JComboBox<String>(new String[] {PercolatorExecutor.V3_01, PercolatorExecutor.V2_10});

	private final JFormattedTextField precursorWindowWidth=new JFormattedTextField(NumberFormat.getNumberInstance());

	private final JComboBox<MassTolerance> precursorTolerance=new JComboBox<MassTolerance>(EncyclopediaParametersPanel.TOLERANCE_VALUES);
	private final JComboBox<MassTolerance> fragmentTolerance=new JComboBox<MassTolerance>(EncyclopediaParametersPanel.TOLERANCE_VALUES);

	private final SpinnerModel minCharge=new SpinnerNumberModel(2, 1, 2, 1);
	private final SpinnerModel maxCharge=new SpinnerNumberModel(3, 2, 4, 1);
	private final SpinnerModel maxMissedCleavage=new SpinnerNumberModel(1, 0, 3, 1);
	private final SpinnerModel numberOfJobs=new SpinnerNumberModel(numberOfCores, 1, numberOfCores, 1);
	private final SpinnerModel numberOfQuantitativeIons=new SpinnerNumberModel(5, 1, 100, 1);
	private final SpinnerModel minNumOfQuantitativeIons=new SpinnerNumberModel(3, 0, 100, 1);
	private final SpinnerModel minQuantitativeIonNumber=new SpinnerNumberModel(3, 0, 100, 1);
	private final JComboBox<String> numberOfExtraDecoyLibraries=new JComboBox<String>(NUMBER_OF_EXTRA_DECOY_ITEMS);

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
		options.add(new LabeledComponent("Data Acquisition Type", acquisition));
		options.add(new LabeledComponent("Precursor Window Width (blank=extract from file)", precursorWindowWidth));
		options.add(new LabeledComponent("Enzyme", enzyme));
		options.add(new LabeledComponent("Fixed", fixed));
		options.add(new LabeledComponent("Fragmentation", fragType));
		options.add(new LabeledComponent("Precursor Mass Tolerance", precursorTolerance));
		options.add(new LabeledComponent("Fragment Mass Tolerance", fragmentTolerance));
		options.add(new LabeledComponent("Maximum Missed Cleavage", new JSpinner(maxMissedCleavage)));
		options.add(new LabeledComponent("Percolator Version", percolatorVersion));
		options.add(new LabeledComponent("Number of Quantitative Ions", new JSpinner(numberOfQuantitativeIons)));
		//options.add(new LabeledComponent("Minimum Number of Quantitative Ions", new JSpinner(minNumOfQuantitativeIons)));
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

	public String getProgramName() {
		return programName;
	}
	
	public String getProgramShortDescription() {
		return programShortDescription;
	}
	
	@Override
	public String getCitation() {
		return "This is a <a href=https://sites.google.com/a/uw.edu/maccoss/>MacCoss Lab</a> project from the University of Washington, <a href=http://www.gs.washington.edu/>Department of Genome Sciences</a>. For more information please contact Brian Searle (searleb@uw.edu).";
	}

	@Override
	public String getAboutMessage() {
		return "<b>Walnut is just like PeCAn, but with more wrinkles and slightly more bitter.";
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
		
		ArrayList<FastaPeptideEntry> targets=null;
		if (targetFile!=null&&!targetFile.equals(fastaFile)) {
			Logger.logLine("Reading targets from ["+targetFile.getName()+"]");
			targets=new ArrayList<FastaPeptideEntry>();
			
			ArrayList<FastaEntryInterface> targetProteins=FastaReader.readFasta(targetFile);
			for (FastaEntryInterface entry : targetProteins) {
				ArrayList<String> peptides=parameters.getEnzyme().digestProtein(entry, parameters.getMinPeptideLength(), parameters.getMaxPeptideLength(), parameters.getMaxMissedCleavages(), parameters.getAAConstants());
				for (String peptide : peptides) {
					FastaPeptideEntry pe=entry.getSubEntry(peptide);
					targets.add(pe);
				}
			}
		}
		
		PercolatorExecutionData percolatorFiles=PecanJobData.getPercolatorExecutionData(diaFile, fastaFile, parameters);
		PecanScoringFactory factory=new PecanOneScoringFactory(parameters, percolatorFiles.getInputTSV());
		return new PecanJob(processor, new PecanJobData(Optional.ofNullable(targets), diaFile, fastaFile, percolatorFiles, factory));
	}

	public PecanSearchParameters getParameters() {
		DataAcquisitionType dataAcquisitionType=DataAcquisitionType.getAcquisitionType((String)acquisition.getSelectedItem());
		DigestionEnzyme digestionEnzyme=DigestionEnzyme.getEnzyme((String)enzyme.getSelectedItem());
		AminoAcidConstants aaConstants=AminoAcidConstants.getConstants((String)fixed.getSelectedItem(), new ModificationMassMap());
		FragmentationType fragmentation=FragmentationType.getFragmentationType((String)fragType.getSelectedItem());
		MassTolerance precursorPPMValue=(MassTolerance)precursorTolerance.getSelectedItem();
		MassTolerance fragmentPPMValue=(MassTolerance)fragmentTolerance.getSelectedItem();
		byte minChargeValue=((Number)minCharge.getValue()).byteValue();
		byte maxChargeValue=((Number)maxCharge.getValue()).byteValue();
		byte maxMissedCleavageValue=((Number)maxMissedCleavage.getValue()).byteValue();
		Number value=(Number)precursorWindowWidth.getValue();
		float precursorWindowWidthValue=value==null?-1.0f:value.floatValue();
		int numberOfJobsValue=((Integer)numberOfJobs.getValue());
		int numberOfQuantitativeIonsValue=((Integer)numberOfQuantitativeIons.getValue());
		int minNumOfQuantitativeIonsValue=((Integer)minNumOfQuantitativeIons.getValue());
		int minQuantitativeIonNumberValue=((Integer)minQuantitativeIonNumber.getValue());
		float numberOfExtraDecoyLibrariesValue=NUMBER_OF_EXTRA_DECOY_VALUES[((Integer)numberOfExtraDecoyLibraries.getSelectedIndex())];
		
		boolean isPercolatorTwo=PercolatorExecutor.V2_10.equals(percolatorVersion.getSelectedItem());
		PecanSearchParameters parameters=new PecanSearchParameters(aaConstants, fragmentation, precursorPPMValue, fragmentPPMValue, digestionEnzyme, isPercolatorTwo?2:3,
				maxMissedCleavageValue, minChargeValue, maxChargeValue, dataAcquisitionType, precursorWindowWidthValue, numberOfJobsValue, numberOfQuantitativeIonsValue, minNumOfQuantitativeIonsValue, minQuantitativeIonNumberValue, numberOfExtraDecoyLibrariesValue, searchPanel.isAlignedBetweenFiles(), true);
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
		
		acquisition.setSelectedItem(DataAcquisitionType.toName(params.getDataAcquisitionType()));
		enzyme.setSelectedItem(params.getEnzyme().getName());
		fixed.setSelectedItem(AminoAcidConstants.toName(params.getAAConstants()));
		fragType.setSelectedItem(FragmentationType.toName(params.getFragType()));
		
		boolean gotIt=false;
		MassTolerance pre=params.getPrecursorTolerance();
		for (int i=0; i<EncyclopediaParametersPanel.TOLERANCE_VALUES.length; i++) {
			if (EncyclopediaParametersPanel.TOLERANCE_VALUES[i].equals(pre)) {
				precursorTolerance.setSelectedIndex(i);
				gotIt=true;
				break;
			}
		}
		if (!gotIt) precursorTolerance.setSelectedIndex(1);
		
		gotIt=false;
		MassTolerance frag=params.getFragmentTolerance();
		for (int i=0; i<EncyclopediaParametersPanel.TOLERANCE_VALUES.length; i++) {
			if (EncyclopediaParametersPanel.TOLERANCE_VALUES[i].equals(frag)) {
				fragmentTolerance.setSelectedIndex(i);
				gotIt=true;
				break;
			}
		}
		if (!gotIt) fragmentTolerance.setSelectedIndex(1);
		
		minCharge.setValue(params.getMinCharge());
		maxCharge.setValue(params.getMaxCharge());
		maxMissedCleavage.setValue(params.getMaxMissedCleavages());
		numberOfJobs.setValue(params.getNumberOfThreadsUsed());
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
		percolatorVersion.setSelectedIndex(params.getPercolatorVersionNumber()==2?1:0);
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
