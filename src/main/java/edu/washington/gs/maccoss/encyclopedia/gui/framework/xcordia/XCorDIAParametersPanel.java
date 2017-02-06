package edu.washington.gs.maccoss.encyclopedia.gui.framework.xcordia;

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

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCordiaSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.ParametersPanelInterface;
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
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassErrorUnitType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

public class XCorDIAParametersPanel extends JPanel implements ParametersPanelInterface {
	private static final long serialVersionUID=1L;
	private static final int numberOfCores=Runtime.getRuntime().availableProcessors();
	private static final String programName="XCorDIA";
	private static final String programShortDescription="XCorDIA Peptide Search";
	public static final ImageIcon smallimage=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/mike_rotate_small_icon.png"));
	public static final ImageIcon image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/mike_rotate_icon.png"));
	public static final String copy="<html><b><p style=\"font-size:16px; font-family: Helvetica, sans-serif\">XCorDIA: Peptide Detection Directly from Data-Independent Acquisition (DIA) MS/MS Data<br></p></b>"
			+ "<p style=\"font-size:10px; font-family: Helvetica, sans-serif\">XCorDIA detects peptides from MZML files, assigns peaks, and calculates various peak features. These features are interpreted by Percolator to identify peptides.";

	private static final String[] NUMBER_OF_EXTRA_DECOY_ITEMS=new String[] {"Normal Target/Decoy", "+10% Extra Decoys", "+20% Extra Decoys", "+50% Extra Decoys", "+100% Extra Decoys (2x Time)"};
	private static final float[] NUMBER_OF_EXTRA_DECOY_VALUES=new float[] {0.0f, 0.1f, 0.2f, 0.5f, 1.0f};
	
	private static final MassTolerance[] TOLERANCE_VALUES=new MassTolerance[] {
			new MassTolerance(5.0, MassErrorUnitType.PPM),  //0
			new MassTolerance(10.0, MassErrorUnitType.PPM), //1
			new MassTolerance(25.0, MassErrorUnitType.PPM), //2
			new MassTolerance(50.0, MassErrorUnitType.PPM), //3
			new MassTolerance(100.0, MassErrorUnitType.PPM),//4
			new MassTolerance(0.4, MassErrorUnitType.AMU),  //5
			new MassTolerance(1.0, MassErrorUnitType.AMU)   //6
	};
	private static final String[] TOLERANCE_NAMES=new String[] {
			TOLERANCE_VALUES[0].toString(), //0
			TOLERANCE_VALUES[1].toString(), //1
			TOLERANCE_VALUES[2].toString(), //2
			TOLERANCE_VALUES[3].toString(), //3
			TOLERANCE_VALUES[4].toString(), //4
			TOLERANCE_VALUES[5].toString(), //5
			TOLERANCE_VALUES[6].toString() //6
	};
	
	private final FileChooserPanel backgroundFasta;
	private final FileChooserPanel targetFasta;
	private final JComboBox<String> acquisition=new JComboBox<String>(new String[] {DataAcquisitionType.toName(DataAcquisitionType.OVERLAPPING_DIA), DataAcquisitionType.toName(DataAcquisitionType.DIA)});
	private final JComboBox<String> enzyme=new JComboBox<String>(new String[] {"Trypsin", "Lys-C", "Lys-N", "Arg-C", "CNBr", "Chymotrypsin", "PepsinA", "No Enzyme"});
	private final JComboBox<String> fixed=new JComboBox<String>(new String[] {"C+57 (Carbamidomethyl)", "C+58 (Carboxymethyl)", "C+46 (MMTS)", "None"});
	private final JComboBox<String> fragType=new JComboBox<String>(new String[] {FragmentationType.toName(FragmentationType.CID), FragmentationType.toName(FragmentationType.YONLY), FragmentationType.toName(FragmentationType.ETD)});
	private final JComboBox<String> precursorTolerance=new JComboBox<String>(TOLERANCE_NAMES);
	private final JComboBox<String> fragmentTolerance=new JComboBox<String>(TOLERANCE_NAMES);

	private final JFormattedTextField precursorWindowWidth=new JFormattedTextField(NumberFormat.getNumberInstance());

	private final SpinnerModel minCharge=new SpinnerNumberModel(2, 1, 2, 1);
	private final SpinnerModel maxCharge=new SpinnerNumberModel(3, 2, 4, 1);
	private final SpinnerModel maxMissedCleavage=new SpinnerNumberModel(1, 0, 3, 1);
	private final SpinnerModel numberOfJobs=new SpinnerNumberModel(numberOfCores, 1, numberOfCores, 1);
	private final SpinnerModel numberOfQuantitativeIons=new SpinnerNumberModel(5, 1, 100, 1);
	private final JComboBox<String> numberOfExtraDecoyLibraries=new JComboBox<String>(NUMBER_OF_EXTRA_DECOY_ITEMS);

	public XCorDIAParametersPanel() {
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
		options.add(new LabeledComponent("Target/Decoy Approach", numberOfExtraDecoyLibraries));
		options.add(new LabeledComponent("Data Acquisition Type", acquisition));
		options.add(new LabeledComponent("Precursor Window Width (blank=extract from file)", precursorWindowWidth));
		options.add(new LabeledComponent("Enzyme", enzyme));
		options.add(new LabeledComponent("Fixed", fixed));
		options.add(new LabeledComponent("Fragmentation", fragType));
		options.add(new LabeledComponent("Precursor (PPM)", precursorTolerance));
		options.add(new LabeledComponent("Fragment (PPM)", fragmentTolerance));
				
		options.add(new LabeledComponent("Maximum Missed Cleavage", new JSpinner(maxMissedCleavage)));
		options.add(new LabeledComponent("Number of Quantitative Ions", new JSpinner(numberOfQuantitativeIons)));
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
	
	@Override
	public void askForSetupFile() {
		backgroundFasta.askForFiles();
	}

	
	public String getProgramName() {
		return programName;
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
		File outputFile=new File(diaFile.getAbsolutePath()+".xcordia.txt");
		File featureFile=new File(outputFile.getAbsolutePath()+".features.txt");
		
		ArrayList<FastaPeptideEntry> targets=null;
		if (targetFile!=null&&!targetFile.equals(fastaFile)) {
			Logger.logLine("Reading targets from ["+targetFile.getName()+"]");
			targets=new ArrayList<FastaPeptideEntry>();
			
			ArrayList<FastaEntryInterface> targetProteins=FastaReader.readFasta(targetFile);
			for (FastaEntryInterface entry : targetProteins) {
				ArrayList<String> peptides=parameters.getEnzyme().digestProtein(entry.getSequence(), parameters.getMinPeptideLength(), parameters.getMaxPeptideLength(), parameters.getMaxMissedCleavages());
				for (String peptide : peptides) {
					FastaPeptideEntry pe=entry.getSubEntry(peptide);
					targets.add(pe);
				}
			}
		}
		
		XCorDIAOneScoringFactory factory=new XCorDIAOneScoringFactory(parameters);
		return new XCorDIAJob(processor, new XCorDIAJobData(Optional.ofNullable(targets), diaFile, fastaFile, featureFile, outputFile, factory));
	}

	public XCordiaSearchParameters getParameters() {
		DataAcquisitionType dataAcquisitionType=DataAcquisitionType.getAcquisitionType((String)acquisition.getSelectedItem());
		DigestionEnzyme digestionEnzyme=DigestionEnzyme.getEnzyme((String)enzyme.getSelectedItem());
		AminoAcidConstants aaConstants=AminoAcidConstants.getConstants((String)fixed.getSelectedItem());
		FragmentationType fragmentation=FragmentationType.getFragmentationType((String)fragType.getSelectedItem());
		MassTolerance precursorPPMValue=TOLERANCE_VALUES[precursorTolerance.getSelectedIndex()];
		MassTolerance fragmentPPMValue=TOLERANCE_VALUES[fragmentTolerance.getSelectedIndex()];
		byte minChargeValue=((Number)minCharge.getValue()).byteValue();
		byte maxChargeValue=((Number)maxCharge.getValue()).byteValue();
		byte maxMissedCleavageValue=((Number)maxMissedCleavage.getValue()).byteValue();
		Number value=(Number)precursorWindowWidth.getValue();
		float precursorWindowWidthValue=value==null?-1.0f:value.floatValue();
		int numberOfJobsValue=((Integer)numberOfJobs.getValue());
		int numberOfQuantitativeIonsValue=((Integer)numberOfQuantitativeIons.getValue());
		float numberOfExtraDecoyLibrariesValue=NUMBER_OF_EXTRA_DECOY_VALUES[((Integer)numberOfExtraDecoyLibraries.getSelectedIndex())];
		XCordiaSearchParameters parameters=new XCordiaSearchParameters(aaConstants, fragmentation, precursorPPMValue, fragmentPPMValue, digestionEnzyme,
				maxMissedCleavageValue, minChargeValue, maxChargeValue, dataAcquisitionType, precursorWindowWidthValue, numberOfJobsValue, numberOfQuantitativeIonsValue, numberOfExtraDecoyLibrariesValue);
		return parameters;
	}
	
	public void setParameters(XCordiaSearchParameters params, String fastaFileName, String targetFileName) {
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
		MassTolerance pre=params.getPrecursorTolerance();
		boolean gotIt=false;
		for (int i=0; i<TOLERANCE_VALUES.length; i++) {
			if (TOLERANCE_VALUES[i].equals(pre)) {
				precursorTolerance.setSelectedIndex(i);
				gotIt=true;
				break;
			}
		}
		if (!gotIt) precursorTolerance.setSelectedIndex(0);
		
		gotIt=false;
		fragmentTolerance.setSelectedIndex(0);
		MassTolerance frag=params.getFragmentTolerance();
		for (int i=0; i<TOLERANCE_VALUES.length; i++) {
			if (TOLERANCE_VALUES[i].equals(frag)) {
				System.out.println("GOT: "+TOLERANCE_VALUES[i]);
				fragmentTolerance.setSelectedIndex(i);
				gotIt=true;
				break;
			}
		}
		if (!gotIt) fragmentTolerance.setSelectedIndex(0);
		
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
