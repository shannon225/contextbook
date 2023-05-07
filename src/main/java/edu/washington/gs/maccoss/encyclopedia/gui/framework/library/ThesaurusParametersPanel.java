package edu.washington.gs.maccoss.encyclopedia.gui.framework.library;

import java.awt.BorderLayout;
import java.awt.Color;
import java.io.File;
import java.text.NumberFormat;
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
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorVersion;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ScoringBreadthType;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ThesaurusJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ThesaurusSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.ParametersPanelInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessorTableModel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LabeledComponent;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.jobs.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.jobs.ThesaurusJob;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassErrorUnitType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class ThesaurusParametersPanel extends JPanel implements ParametersPanelInterface {
	
	private static final long serialVersionUID=1L;
	private static final int numberOfCores=Runtime.getRuntime().availableProcessors();
	
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
	
	public static final ScoringBreadthType[] CASiL_SEARCH_TYPES=new ScoringBreadthType[] {
			ScoringBreadthType.ENTIRE_RT_WINDOW,
			ScoringBreadthType.RECALIBRATED_20_PERCENT,
			ScoringBreadthType.RECALIBRATED_PEAK_WIDTH,
			ScoringBreadthType.UNCALIBRATED_20_PERCENT,
			ScoringBreadthType.UNCALIBRATED_PEAK_WIDTH
	};
	
	private static final ImageIcon smallimage=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/thesaurus_small_icon.png"));
	private static final ImageIcon image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/thesaurus_icon.png"));
	private static final String programName="Thesaurus";
	private static final String programShortDescription="Phosphopeptide Positional Isomer Search Engine";
	private static final String copy="<html><b><p style=\"font-size:16px; font-family: Helvetica, sans-serif\">Thesaurus: Software for Detecting Positional Phosphopeptide Isomers from Data-Independent Acquisition (DIA) MS/MS Data<br></p></b>"
			+ "<p style=\"font-size:10px; font-family: Helvetica, sans-serif\">Thesaurus extracts peptide fragmentation chromatograms from MZML files, matches them to spectra in libraries, and calculates various scoring features. Matches are localized and alternate positional isomers are explored. These isomers are interpreted by Percolator to identify site-specific peptides.";

	private final FileChooserPanel backgroundFasta=new FileChooserPanel(null, "Background", new SimpleFilenameFilter(".fas", ".fasta"), true);
	private final FileChooserPanel libraryFileChooser;
	private final JComboBox<String> enzyme=new JComboBox<String>(new String[] {"Trypsin", "Lys-C", "Lys-N", "Arg-C", "CNBr", "Chymotrypsin", "Pepsin A"});
	private final JComboBox<String> fragType=new JComboBox<String>(new String[] {FragmentationType.toName(FragmentationType.CID), FragmentationType.toName(FragmentationType.HCD), //FragmentationType.toName(FragmentationType.SILAC), 
			FragmentationType.toName(FragmentationType.ETD)});
	private final JComboBox<PercolatorVersion> percolatorVersion=new JComboBox<PercolatorVersion>(PercolatorVersion.VALID_VERSIONS);
	
	private final JFormattedTextField precursorWindowWidth=new JFormattedTextField(NumberFormat.getNumberInstance()); // not displayed anymore

	private final JComboBox<MassTolerance> precursorTolerance=new JComboBox<MassTolerance>(TOLERANCE_VALUES);
	private final JComboBox<MassTolerance> fragmentTolerance=new JComboBox<MassTolerance>(TOLERANCE_VALUES);
	private final JComboBox<MassTolerance> libraryTolerance=new JComboBox<MassTolerance>(TOLERANCE_VALUES);
	private final SpinnerModel numberOfJobs=new SpinnerNumberModel(numberOfCores, 1, numberOfCores, 1);
	private final JComboBox<ScoringBreadthType> searchBreadthType=new JComboBox<>(CASiL_SEARCH_TYPES);
	private final JComboBox<PeptideModification> modificationType=new JComboBox<>(AminoAcidConstants.getDefaultLocalizationModifications().toArray(new PeptideModification[0]));
	private final SpinnerModel numberOfQuantitativeIons=new SpinnerNumberModel(5, 1, 100, 1);
	private final SpinnerModel percolatorThreshold=new SpinnerNumberModel(0.01, 0.001, 0.1, 0.001);
	private final SpinnerModel minNumOfQuantitativeIons=new SpinnerNumberModel(3, 0, 100, 1);
	private final JTextField additionalCommandLineOptions=new JTextField();
	private final SearchPanel searchPanel;

	public ThesaurusParametersPanel(SearchPanel searchPanel) {
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
		options.add(new LabeledComponent("Modification Type", modificationType));
		options.add(new LabeledComponent("Localization Strategy", searchBreadthType));
		options.add(new LabeledComponent("Enzyme", enzyme));
		options.add(new LabeledComponent("Fragmentation", fragType));
		options.add(new LabeledComponent("Precursor Mass Tolerance", precursorTolerance));
		options.add(new LabeledComponent("Fragment Mass Tolerance", fragmentTolerance));
		options.add(new LabeledComponent("Library Mass Tolerance", libraryTolerance));
		options.add(new LabeledComponent("Percolator Version", percolatorVersion));
		options.add(new LabeledComponent("Percolator FDR threshold", new JSpinner(percolatorThreshold)));
		//options.add(new LabeledComponent("Number of Quantitative Ions", new JSpinner(numberOfQuantitativeIons)));
		options.add(new LabeledComponent("Minimum Number of Well-Shaped Ions", new JSpinner(minNumOfQuantitativeIons)));
		options.add(new LabeledComponent("Number of Cores", new JSpinner(numberOfJobs)));
		options.add(new LabeledComponent("Additonal Command Line Options", additionalCommandLineOptions));

		this.add(options, BorderLayout.CENTER);
	}

	public ProgramType getProgram() {
		return ProgramType.CASiL;
	}
	
	public String getProgramShortDescription() {
		return programShortDescription;
	}
	
	public ImageIcon getSmallImage() {
		return smallimage;
	}

	@Override
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
	
	@Override
	public boolean canCreateChromatogramLibraries() {
		return true;
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
	public void getJob(File diaFile, JobProcessorTableModel model) {
		SearchParameters parameters=getParameters();
		File libraryFile=libraryFileChooser.getFile();
		File fastaFile=getBackgroundFastaFile();
		if (libraryFile==null) return;
		SearchJob job=getJob(diaFile, libraryFile, fastaFile, parameters);

		if (job!=null) {
			model.addJob(job);
		}
	}

	private static HashMap<File, LibraryInterface> libraries=new HashMap<File, LibraryInterface>();
	static SearchJob getJob(File diaFile, File libraryFile, File fastaFile, SearchParameters parameters) {
		File outputFile=new File(diaFile.getAbsolutePath()+ThesaurusJobData.OUTPUT_FILE_SUFFIX);
		
		LibraryInterface library=libraries.get(libraryFile);
		if (library==null) {
			library=BlibToLibraryConverter.getFile(libraryFile, fastaFile, parameters);
			libraries.put(libraryFile, library);
		}
		
		// Thesaurus only tested with V1
		LibraryScoringFactory factory=new EncyclopediaOneScoringFactory(parameters);
		ThesaurusJobData job=new ThesaurusJobData(diaFile, library, outputFile, fastaFile, factory);
		return new ThesaurusJob(job);
	}

	public ThesaurusSearchParameters getParameters() {
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
		float targetWindowCenter=-1f;
		int numberOfQuantitativeIonsValue=((Integer)numberOfQuantitativeIons.getValue());
		int minNumOfQuantitativeIonsValue=((Integer)minNumOfQuantitativeIons.getValue());
		ScoringBreadthType CASiLSearchBreadthType=(ScoringBreadthType)searchBreadthType.getSelectedItem();
		PeptideModification modification=(PeptideModification)modificationType.getSelectedItem();
		float percolatorThresholdValue=((Number)percolatorThreshold.getValue()).floatValue();
		
		boolean considerRearrangement=false;
		
		ThesaurusSearchParameters parameters=new ThesaurusSearchParameters(
				aaConstants,
				fragmentation,
				precursorValue,
				0.0,
				0.0,
				fragmentValue,
				0.0,
				libraryFragmentValue,
				digestionEnzyme,
				percolatorThresholdValue,
				percolatorThresholdValue,
				percolator,
				PercolatorExecutor.DEFAULT_TRAINING_SET_SIZE,
				PercolatorExecutor.DEFAULT_TRAINING_THRESHOLD,
				PercolatorExecutor.DEFAULT_TRAINING_ITERATIONS,
				dataAcquisitionType,
				numberOfJobsValue,
				25f,
				targetWindowCenter,
				precursorWindowWidthValue,
				numberOfQuantitativeIonsValue,
				minNumOfQuantitativeIonsValue,
				-1,
				0.0f,
				false,
				modification,
				CASiLSearchBreadthType,
				0.0f,
				true,
				false,
				false,
				false,
				considerRearrangement,
				Optional.empty(),
				Optional.empty(),
				true,
				true,
				false,
				false,
				false
		);

		String cmds=additionalCommandLineOptions.getText();
		HashMap<String, String> params=parameters.toParameterMap();
		params.putAll(CommandLineParser.parseArguments(cmds.split(" ")));
		parameters=ThesaurusSearchParameters.parseParameters(params);
		
		return parameters;
	}
	
	public void setParameters(ThesaurusSearchParameters params, String libraryFileName, String fastaFileName) {
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
		
		searchBreadthType.setSelectedItem(params.getScoringBreadthType());
		if (params.getLocalizingModification().isPresent()) {
			modificationType.setSelectedItem(params.getLocalizingModification().get());
		}
		
		numberOfJobs.setValue(params.getNumberOfThreadsUsed());
		if (params.getPrecursorWindowSize()>0) {
			precursorWindowWidth.setValue(params.getPrecursorWindowSize());
		} else {
			precursorWindowWidth.setValue(-1);
		}
		percolatorVersion.setSelectedItem(params.getPercolatorVersionNumber());
		numberOfQuantitativeIons.setValue(params.getNumberOfQuantitativePeaks());
		minNumOfQuantitativeIons.setValue(params.getMinNumOfQuantitativePeaks());
		percolatorThreshold.setValue(new Double(params.getPercolatorThreshold()));
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
