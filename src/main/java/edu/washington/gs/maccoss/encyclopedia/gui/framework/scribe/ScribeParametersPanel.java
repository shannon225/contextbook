package edu.washington.gs.maccoss.encyclopedia.gui.framework.scribe;

import java.awt.BorderLayout;
import java.awt.Color;
import java.io.File;
import java.text.NumberFormat;
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
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorVersion;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ScoringBreadthType;
import edu.washington.gs.maccoss.encyclopedia.algorithms.scribe.ScribeJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.scribe.ScribeScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.scribe.ScribeSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.ParametersPanelInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessorTableModel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LabeledComponent;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingJob;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassErrorUnitType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class ScribeParametersPanel extends JPanel implements ParametersPanelInterface {
	
	private static final long serialVersionUID=1L;
	private static final int numberOfCores=Runtime.getRuntime().availableProcessors();
	private static final String[] NUMBER_OF_EXTRA_DECOY_ITEMS=new String[] {"Normal Target/Decoy", "+10% Extra Decoys", "+20% Extra Decoys", "+50% Extra Decoys", "+100% Extra Decoys (2x Time)"};
	private static final float[] NUMBER_OF_EXTRA_DECOY_VALUES=new float[] {0.0f, 0.1f, 0.2f, 0.5f, 1.0f};
	
	public static final MassTolerance[] TOLERANCE_VALUES=new MassTolerance[] {
			new MassTolerance(5.0, MassErrorUnitType.PPM),  //0
			new MassTolerance(10.0, MassErrorUnitType.PPM), //1
			new MassTolerance(25.0, MassErrorUnitType.PPM), //2
			new MassTolerance(50.0, MassErrorUnitType.PPM), //3
			new MassTolerance(100.0, MassErrorUnitType.PPM),//4
			new MassTolerance(0.4, MassErrorUnitType.AMU),  //5
			new MassTolerance(1.0, MassErrorUnitType.AMU),   //6
			new MassTolerance(15000.0, MassErrorUnitType.RESOLUTION),//7
			new MassTolerance(17500.0, MassErrorUnitType.RESOLUTION), //8
			new MassTolerance(30000.0, MassErrorUnitType.RESOLUTION),//9
			new MassTolerance(35000.0, MassErrorUnitType.RESOLUTION), //10
			new MassTolerance(60000.0, MassErrorUnitType.RESOLUTION),//11
	};
	
	private static final ImageIcon smallimage=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/scribe_small_icon.png"));
	private static final ImageIcon image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/scribe_icon.png"));
	private static final String programName="Scribe";
	private static final String programShortDescription="Scribe Library Search";
	private static final String copy="<html><b><p style=\"font-size:16px; font-family: Helvetica, sans-serif\">Scribe: Library Searching from Data-Dependent Acquisition (DDA) MS/MS Data<br></p></b>"
			+ "<p style=\"font-size:10px; font-family: Helvetica, sans-serif\">Scribe extracts peptide fragmentation spectra from MZML files, matches them to spectra in libraries, and calculates various scoring features. These features are interpreted by Percolator to identify peptides.";

	private final FileChooserPanel backgroundFasta=new FileChooserPanel(null, "Background", new SimpleFilenameFilter(".fas", ".fasta"), true);
	private final FileChooserPanel libraryFileChooser;
	private final JComboBox<String> enzyme=new JComboBox<String>(new String[] {"Trypsin", "Glu-C", "Lys-C", "Arg-C", "Asp-N", "Lys-N", "CNBr", "Chymotrypsin", "Pepsin A", "No Enzyme"});
	private final JComboBox<String> fragType=new JComboBox<String>(new String[] {FragmentationType.toName(FragmentationType.CID), FragmentationType.toName(FragmentationType.HCD), FragmentationType.toName(FragmentationType.ETD)});
	private final JComboBox<PercolatorVersion> percolatorVersion=new JComboBox<PercolatorVersion>(PercolatorVersion.VALID_VERSIONS);

	private final JFormattedTextField precursorWindowWidth=new JFormattedTextField(NumberFormat.getNumberInstance()); // not displayed anymore

	private final JComboBox<MassTolerance> precursorTolerance=new JComboBox<MassTolerance>(TOLERANCE_VALUES);
	private final JComboBox<MassTolerance> fragmentTolerance=new JComboBox<MassTolerance>(TOLERANCE_VALUES);
	private final JComboBox<MassTolerance> libraryTolerance=new JComboBox<MassTolerance>(TOLERANCE_VALUES);
	private final SpinnerModel numberOfJobs=new SpinnerNumberModel(numberOfCores, 1, numberOfCores, 1);
	private final SpinnerModel numberOfQuantitativeIons=new SpinnerNumberModel(5, 1, 100, 1);
	private final SpinnerModel minNumOfQuantitativeIons=new SpinnerNumberModel(3, 0, 100, 1);
	private final JComboBox<String> numberOfExtraDecoyLibraries=new JComboBox<String>(NUMBER_OF_EXTRA_DECOY_ITEMS);
	private final JTextField additionalCommandLineOptions=new JTextField();

	private final SearchPanel searchPanel;
	public ScribeParametersPanel(SearchPanel searchPanel) {
		super(new BorderLayout());
		this.searchPanel=searchPanel;

		JPanel top=new JPanel(new BorderLayout());
		top.add(new JLabel(getImage()), BorderLayout.WEST);
		JEditorPane editor=new JEditorPane("text/html", getCopy());
		editor.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
		top.add(editor, BorderLayout.CENTER);
		top.setBackground(Color.white);
		this.add(top, BorderLayout.NORTH);
		
		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(new LabeledComponent("<p style=\"font-size:12px; font-family: Helvetica, sans-serif\"><b>Parameters", new JLabel()));

		libraryFileChooser=new FileChooserPanel(null, "Library", new SimpleFilenameFilter(LibraryFile.DLIB, LibraryFile.ELIB), true);
		options.add(libraryFileChooser);
		options.add(backgroundFasta);
		options.add(new LabeledComponent("Target/Decoy Approach", numberOfExtraDecoyLibraries));
		options.add(new LabeledComponent("Enzyme", enzyme));
		options.add(new LabeledComponent("Fragmentation", fragType));
		options.add(new LabeledComponent("Precursor Mass Tolerance", precursorTolerance));
		options.add(new LabeledComponent("Fragment Mass Tolerance", fragmentTolerance));
		options.add(new LabeledComponent("Library Mass Tolerance", libraryTolerance));
		options.add(new LabeledComponent("Percolator Version", percolatorVersion));
		options.add(new LabeledComponent("Number of Quantitative Ions", new JSpinner(numberOfQuantitativeIons)));
		options.add(new LabeledComponent("Minimum Number of Quantitative Ions", new JSpinner(minNumOfQuantitativeIons)));
		options.add(new LabeledComponent("Number of Cores", new JSpinner(numberOfJobs)));
		options.add(new LabeledComponent("Additonal Command Line Options", additionalCommandLineOptions));

		this.add(options, BorderLayout.CENTER);
	}
	
	public String getProgramName() {
		return programName;
	}
	
	public String getProgramShortDescription() {
		return programShortDescription;
	}

	public ProgramType getProgram() {
		return ProgramType.Scribe;
	}
	
	public ImageIcon getSmallImage() {
		return smallimage;
	}
	
	public ImageIcon getImage() {
		return image;
	}
	
	public String getCopy() {
		return copy;
	}
	
	@Override
	public void askForSetupFile() {
		libraryFileChooser.askForFiles();
	}

	@Override
	public File getBackgroundFastaFile() {
		return backgroundFasta.getFile();
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.gui.pecan.ParametersPanelInterface#canLoadData()
	 */
	@Override
	public Optional<String> canLoadData() {
		if (libraryFileChooser.getFile()==null) {
			return Optional.of("Please load a library file first!");	
		}
		if (getBackgroundFastaFile()==null) {
			return Optional.of("Please load a fasta file first!");	
		}
		return Optional.empty();
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.gui.pecan.ParametersPanelInterface#getJob(java.io.File, edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessorTableModel)
	 */
	@Override
	public SwingJob getJob(File diaFile, JobProcessorTableModel model) {
		SearchParameters parameters=getParameters();
		File libraryFile=libraryFileChooser.getFile();
		File fastaFile=getBackgroundFastaFile();
		if (libraryFile==null) return null;
		SearchJob job=getJob(diaFile, fastaFile, libraryFile, model, parameters);

		if (job!=null) {
			model.addJob(job);
		}
		return job;
	}

	private static HashMap<File, LibraryInterface> libraries=new HashMap<File, LibraryInterface>();
	static SearchJob getJob(File diaFile, File fastaFile, File libraryFile, JobProcessor processor, SearchParameters parameters) {
		
		LibraryInterface library=libraries.get(libraryFile);
		if (library==null) {
			library=BlibToLibraryConverter.getFile(libraryFile, fastaFile, parameters);
			libraries.put(libraryFile, library);
		}
		
		ScribeScoringFactory factory=new ScribeScoringFactory(parameters);
		ScribeJobData job=new ScribeJobData(diaFile, fastaFile, library, factory);
		return new ScribeJob(processor, job);
	}

	public ScribeSearchParameters getParameters() {
		DataAcquisitionType dataAcquisitionType=DataAcquisitionType.DIA;
		DigestionEnzyme digestionEnzyme=DigestionEnzyme.getEnzyme((String)enzyme.getSelectedItem());
		AminoAcidConstants aaConstants=new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());
		FragmentationType fragmentation=FragmentationType.getFragmentationType((String)fragType.getSelectedItem());
		MassTolerance precursorValue=(MassTolerance)precursorTolerance.getSelectedItem();
		MassTolerance fragmentValue=(MassTolerance)fragmentTolerance.getSelectedItem();
		MassTolerance libraryFragmentValue=(MassTolerance)libraryTolerance.getSelectedItem();
		int numberOfJobsValue=((Integer)numberOfJobs.getValue());
		Number value=(Number)precursorWindowWidth.getValue();
		float precursorWindowWidthValue=value==null?-1.0f:value.floatValue();
		PercolatorVersion percolator=(PercolatorVersion)percolatorVersion.getSelectedItem();
		float numberOfExtraDecoyLibrariesValue=NUMBER_OF_EXTRA_DECOY_VALUES[((Integer)numberOfExtraDecoyLibraries.getSelectedIndex())];
		float targetWindowCenter=-1f;
		int numberOfQuantitativeIonsValue=((Integer)numberOfQuantitativeIons.getValue());
		int minNumOfQuantitativeIonsValue=((Integer)minNumOfQuantitativeIons.getValue());
		float minIntensity=-1.0f;
		Optional<PeptideModification> modificationType=Optional.empty();
		ScribeSearchParameters parameters=new ScribeSearchParameters(
				aaConstants,
				fragmentation,
				precursorValue,
				0.0,
				0.0,
				fragmentValue,
				0.0,
				libraryFragmentValue,
				digestionEnzyme,
				0.01f,
				0.01f,
				percolator,
				PercolatorExecutor.DEFAULT_TRAINING_SET_SIZE,
				PercolatorExecutor.DEFAULT_TRAINING_THRESHOLD,
				dataAcquisitionType,
				numberOfJobsValue,
				25f,
				targetWindowCenter,
				precursorWindowWidthValue,
				numberOfQuantitativeIonsValue,
				minNumOfQuantitativeIonsValue,
				-1,
				minIntensity,
				modificationType,
				ScoringBreadthType.ENTIRE_RT_WINDOW,
				numberOfExtraDecoyLibrariesValue,
				true,
				true,
				false,
				false,
				Optional.empty(),
				Optional.empty(),
				false
		);

		String cmds=additionalCommandLineOptions.getText();
		HashMap<String, String> params=parameters.toParameterMap();
		params.putAll(CommandLineParser.parseArguments(cmds.split(" ")));
		parameters=ScribeSearchParameters.convertFromEncyclopeDIA(SearchParameterParser.parseParameters(params));
		
		return parameters;
	}
	
	public void setParameters(SearchParameters params, String libraryFileName, String fastaFileName) {
		if (libraryFileName!=null) {
			File libraryFile=new File(libraryFileName);
			if (libraryFile.exists()) libraryFileChooser.update(libraryFile);
		}
		if (fastaFileName!=null) {
			File fastaFile=new File(fastaFileName);
			if (fastaFile.exists()) backgroundFasta.update(fastaFile);
		}
		
		enzyme.setSelectedItem(params.getEnzyme().getName());
		fragType.setSelectedItem(FragmentationType.toName(params.getFragType()));
		
		boolean gotIt=false;
		MassTolerance pre=params.getPrecursorTolerance();
		for (int i=0; i<TOLERANCE_VALUES.length; i++) {
			if (TOLERANCE_VALUES[i].equals(pre)) {
				precursorTolerance.setSelectedIndex(i);
				gotIt=true;
				break;
			}
		}
		if (!gotIt) precursorTolerance.setSelectedIndex(1);
		
		gotIt=false;
		MassTolerance frag=params.getFragmentTolerance();
		for (int i=0; i<TOLERANCE_VALUES.length; i++) {
			if (TOLERANCE_VALUES[i].equals(frag)) {
				fragmentTolerance.setSelectedIndex(i);
				gotIt=true;
				break;
			}
		}
		if (!gotIt) fragmentTolerance.setSelectedIndex(1);
		
		gotIt=false;
		MassTolerance lib=params.getLibraryFragmentTolerance();
		for (int i=0; i<TOLERANCE_VALUES.length; i++) {
			if (TOLERANCE_VALUES[i].equals(lib)) {
				libraryTolerance.setSelectedIndex(i);
				gotIt=true;
				break;
			}
		}
		if (!gotIt) libraryTolerance.setSelectedIndex(1);
		
		numberOfJobs.setValue(params.getNumberOfThreadsUsed());
		if (params.getPrecursorWindowSize()>0) {
			precursorWindowWidth.setValue(params.getPrecursorWindowSize());
		} else {
			precursorWindowWidth.setValue(-1);
		}
		percolatorVersion.setSelectedItem(params.getPercolatorVersionNumber());
		int index=Arrays.binarySearch(NUMBER_OF_EXTRA_DECOY_VALUES, params.getNumberOfExtraDecoyLibrariesSearched());
		if (index>=0) {
			numberOfExtraDecoyLibraries.setSelectedIndex(index);
		}
		numberOfQuantitativeIons.setValue(params.getNumberOfQuantitativePeaks());
		minNumOfQuantitativeIons.setValue(params.getMinNumOfQuantitativePeaks());
	}
	
	@Override
	public void savePreferences() {
		try {
			getParameters().savePreferences(libraryFileChooser.getFile(), backgroundFasta.getFile());
		} catch (Exception e) {
			Logger.errorLine("Error writing parameters to disk!");
			Logger.errorException(e);
		}
	}
}
