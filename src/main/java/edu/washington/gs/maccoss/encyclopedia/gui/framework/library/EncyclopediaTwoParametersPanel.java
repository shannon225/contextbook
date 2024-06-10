package edu.washington.gs.maccoss.encyclopedia.gui.framework.library;

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
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaTwoAlignmentLibraryFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaTwoJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorVersion;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ScoringBreadthType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.ParametersPanelInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessorTableModel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LabeledComponent;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.jobs.EncyclopediaTwoJob;
import edu.washington.gs.maccoss.encyclopedia.jobs.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassErrorUnitType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class EncyclopediaTwoParametersPanel extends JPanel implements ParametersPanelInterface {
	
	private static final long serialVersionUID=1L;
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
	
	private static final ImageIcon smallimage=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/encyclopedia_small_icon.png"));
	private static final ImageIcon image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/encyclopedia_icon.png"));
	private static final String programName="EncyclopeDIA";
	private static final String programShortDescription="EncyclopeDIA Library Search";
	private static final String copy="<html><b><p style=\"font-size:16px; font-family: Helvetica, sans-serif\">EncyclopeDIA: Library Searching Directly from Data-Independent Acquisition (DIA) MS/MS Data<br></p></b>"
			+ "<p style=\"font-size:10px; font-family: Helvetica, sans-serif\">EncyclopeDIA extracts peptide fragmentation chromatograms from MZML files, matches them to spectra in libraries, and calculates various scoring features. These features are interpreted by Percolator to identify peptides.";

	private final FileChooserPanel backgroundFasta=new FileChooserPanel(null, "Background", new SimpleFilenameFilter(".fas", ".fasta"), true);
	private final FileChooserPanel prealignmentLibraryFileChooser;
	private final FileChooserPanel libraryFileChooser;
	private final JComboBox<String> fragType=new JComboBox<String>(new String[] {FragmentationType.toName(FragmentationType.CID), FragmentationType.toName(FragmentationType.HCD), //FragmentationType.toName(FragmentationType.SILAC), 
			FragmentationType.toName(FragmentationType.ETD)});
	private final JComboBox<PercolatorVersion> percolatorVersion=new JComboBox<PercolatorVersion>(PercolatorVersion.VALID_VERSIONS);

	private final JFormattedTextField precursorWindowWidth=new JFormattedTextField(NumberFormat.getNumberInstance()); // not displayed anymore

	private final JComboBox<MassTolerance> libraryTolerance=new JComboBox<MassTolerance>(TOLERANCE_VALUES);
	private final SpinnerModel numberOfQuantitativeIons=new SpinnerNumberModel(5, 1, 100, 1);
	private final SpinnerModel minNumOfQuantitativeIons=new SpinnerNumberModel(3, 0, 100, 1);
	private final JComboBox<String> numberOfExtraDecoyLibraries=new JComboBox<String>(NUMBER_OF_EXTRA_DECOY_ITEMS);
	private final JTextField additionalCommandLineOptions=new JTextField();

	private final SearchPanel searchPanel;
	public EncyclopediaTwoParametersPanel(SearchPanel searchPanel) {
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

		prealignmentLibraryFileChooser=new FileChooserPanel(null, "Prealignment Library", new SimpleFilenameFilter(LibraryFile.DLIB, LibraryFile.ELIB), false);
		options.add(prealignmentLibraryFileChooser);
		
		libraryFileChooser=new FileChooserPanel(null, "Library", new SimpleFilenameFilter(LibraryFile.DLIB, LibraryFile.ELIB), true);
		options.add(libraryFileChooser);
		options.add(backgroundFasta);
		options.add(new LabeledComponent("Target/Decoy Approach", numberOfExtraDecoyLibraries));
		options.add(new LabeledComponent("Library Mass Tolerance", libraryTolerance));
		options.add(new LabeledComponent("Number of Quantitative Ions", new JSpinner(numberOfQuantitativeIons)));
		options.add(new LabeledComponent("Minimum Number of Quantitative Ions", new JSpinner(minNumOfQuantitativeIons)));
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
		return ProgramType.EncyclopeDIA;
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
	public boolean canCreateChromatogramLibraries() {
		return true;
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
	public void getJob(File diaFile, JobProcessorTableModel model) {
		SearchParameters parameters=getParameters();
		File prealignmentLibraryFile=prealignmentLibraryFileChooser.getFile();
		File libraryFile=libraryFileChooser.getFile();
		File fastaFile=getBackgroundFastaFile();
		if (prealignmentLibraryFile==null) {
			Optional<File> optFile=EncyclopediaTwoAlignmentLibraryFactory.getPreAlignmentFile(searchPanel.getEnzyme());
			if (!optFile.isPresent()) {
				Logger.errorLine("Sorry, no pre-alignment library for that enzyme. You must select a library or change enzyme!");
				return;
			}
			prealignmentLibraryFile=optFile.get();
			Logger.logLine("Using prealignment library: "+prealignmentLibraryFile.getName());
		}
		if (libraryFile==null) return;
		SearchJob job=getJob(diaFile, fastaFile, prealignmentLibraryFile, libraryFile, parameters);

		if (job!=null) {
			model.addJob(job);
		}
	}

	private static HashMap<File, LibraryInterface> libraries=new HashMap<File, LibraryInterface>();
	static SearchJob getJob(File diaFile, File fastaFile, File prealignmentLibraryFile, File libraryFile, SearchParameters parameters) {
		
		LibraryInterface prelib=libraries.get(prealignmentLibraryFile);
		if (prelib==null) {
			prelib=BlibToLibraryConverter.getFile(prealignmentLibraryFile, fastaFile, parameters);
			libraries.put(prealignmentLibraryFile, prelib);
		}
		LibraryInterface library=libraries.get(libraryFile);
		if (library==null) {
			library=BlibToLibraryConverter.getFile(libraryFile, fastaFile, parameters);
			libraries.put(libraryFile, library);
		}

		LibraryScoringFactory factory=EncyclopediaScoringFactory.getDefaultScoringFactory(parameters);
		EncyclopediaTwoJobData job=new EncyclopediaTwoJobData(diaFile, fastaFile, prelib, library, factory);
		return new EncyclopediaTwoJob(job);
	}

	public SearchParameters getParameters() {
		DataAcquisitionType dataAcquisitionType=DataAcquisitionType.DIA;
		DigestionEnzyme digestionEnzyme=searchPanel.getEnzyme();
		AminoAcidConstants aaConstants=new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());
		FragmentationType fragmentation=FragmentationType.getFragmentationType((String)fragType.getSelectedItem());
		MassTolerance precursorValue=searchPanel.getInstrument().getPrecursorTolerance();
		MassTolerance fragmentValue=searchPanel.getInstrument().getFragmentTolerance();
		MassTolerance libraryFragmentValue=(MassTolerance)libraryTolerance.getSelectedItem();
		int numberOfJobsValue=searchPanel.getNumberOfJobs();
		Number value=(Number)precursorWindowWidth.getValue();
		float precursorWindowWidthValue=value==null?-1.0f:value.floatValue();
		PercolatorVersion percolator=(PercolatorVersion)percolatorVersion.getSelectedItem();
		float numberOfExtraDecoyLibrariesValue=NUMBER_OF_EXTRA_DECOY_VALUES[((Integer)numberOfExtraDecoyLibraries.getSelectedIndex())];
		float targetWindowCenter=-1f;
		int numberOfQuantitativeIonsValue=((Integer)numberOfQuantitativeIons.getValue());
		int minNumOfQuantitativeIonsValue=((Integer)minNumOfQuantitativeIons.getValue());
		float minIntensity=-1.0f;
		float IITNumberOfIonsThreshold=-1.0f;
		Optional<PeptideModification> modificationType=Optional.empty();
		SearchParameters parameters=new SearchParameters(
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
				true,
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
				SearchParameters.DEFAULT_MAX_NUM_FRAGMENT_IONS,
				minIntensity,
				IITNumberOfIonsThreshold,
				modificationType,
				ScoringBreadthType.ENTIRE_RT_WINDOW,
				numberOfExtraDecoyLibrariesValue,
				true,
				true,
				-1.0f,
				SearchParameters.DEFAULT_MIN_NUM_INTEGRATED_RT_POINTS,
				false,
				false,
				Optional.empty(),
				Optional.empty(),
				true,
				true,
				false,
				false,
				false,
				searchPanel.getInstrument(),
				false
		);

		String cmds=additionalCommandLineOptions.getText();
		HashMap<String, String> params=parameters.toParameterMap();
		params=searchPanel.getInstrument().overwriteParameters(params);
		params.putAll(CommandLineParser.parseArguments(cmds.split(" ")));
		parameters=SearchParameterParser.parseParameters(params);
		
		return parameters;
	}
	
	public void setParameters(SearchParameters params, String preLibraryFileName, String libraryFileName, String fastaFileName) {
		if (preLibraryFileName!=null) {
			File preLibraryFile=new File(preLibraryFileName);
			if (preLibraryFile.exists()) prealignmentLibraryFileChooser.update(preLibraryFile);
		}
		if (libraryFileName!=null) {
			File libraryFile=new File(libraryFileName);
			if (libraryFile.exists()) libraryFileChooser.update(libraryFile);
		}
		if (fastaFileName!=null) {
			File fastaFile=new File(fastaFileName);
			if (fastaFile.exists()) backgroundFasta.update(fastaFile);
		}
		
		boolean gotIt=false;
		MassTolerance lib=params.getLibraryFragmentTolerance();
		for (int i=0; i<TOLERANCE_VALUES.length; i++) {
			if (TOLERANCE_VALUES[i].equals(lib)) {
				libraryTolerance.setSelectedIndex(i);
				gotIt=true;
				break;
			}
		}
		if (!gotIt) libraryTolerance.setSelectedIndex(1);
		
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
