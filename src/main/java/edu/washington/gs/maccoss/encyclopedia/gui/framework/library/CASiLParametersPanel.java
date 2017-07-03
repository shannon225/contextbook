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
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.CASiLJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.CASiLSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ScoringBreadthType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
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
import gnu.trove.map.hash.TCharDoubleHashMap;

public class CASiLParametersPanel extends JPanel implements ParametersPanelInterface {
	
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
	private static final String programShortDescription="Chromatogram-Aligned Site Localizing Search Engine";
	private static final String copy="<html><b><p style=\"font-size:16px; font-family: Helvetica, sans-serif\">Thesaurus: Chromatogram-Aligned Site Localizing Search Engine for Data-Independent Acquisition (DIA) MS/MS Data<br></p></b>"
			+ "<p style=\"font-size:10px; font-family: Helvetica, sans-serif\">Thesaurus extracts peptide fragmentation chromatograms from MZML files, matches them to spectra in libraries, and calculates various scoring features. Matches are localized and alternate positional isomers are explored. These isomers are interpreted by Percolator to identify site-specific peptides at an estimated 1% FDR threshold.";
	
	private final FileChooserPanel libraryFileChooser;
	private final JComboBox<String> acquisition=new JComboBox<String>(new String[] {DataAcquisitionType.toName(DataAcquisitionType.OVERLAPPING_DIA), DataAcquisitionType.toName(DataAcquisitionType.DIA), DataAcquisitionType.toName(DataAcquisitionType.DDA)});
	private final JComboBox<String> enzyme=new JComboBox<String>(new String[] {"Trypsin", "Lys-C", "Lys-N", "Arg-C", "CNBr", "Chymotrypsin", "Pepsin A"});
	private final JComboBox<String> fragType=new JComboBox<String>(new String[] {FragmentationType.toName(FragmentationType.CID), FragmentationType.toName(FragmentationType.YONLY), FragmentationType.toName(FragmentationType.ETD)});
	private final JComboBox<String> percolatorVersion=new JComboBox<String>(new String[] {PercolatorExecutor.V3_01, PercolatorExecutor.V2_10});

	private final JFormattedTextField precursorWindowWidth=new JFormattedTextField(NumberFormat.getNumberInstance()); // not displayed anymore

	private final JComboBox<MassTolerance> precursorTolerance=new JComboBox<MassTolerance>(TOLERANCE_VALUES);
	private final JComboBox<MassTolerance> fragmentTolerance=new JComboBox<MassTolerance>(TOLERANCE_VALUES);
	private final JComboBox<MassTolerance> libraryTolerance=new JComboBox<MassTolerance>(TOLERANCE_VALUES);
	private final SpinnerModel numberOfJobs=new SpinnerNumberModel(numberOfCores, 1, numberOfCores, 1);
	private final JComboBox<ScoringBreadthType> searchBreadthType=new JComboBox<>(CASiL_SEARCH_TYPES);
	private final JComboBox<PeptideModification> modificationType=new JComboBox<>(PeptideModification.MODIFICATIONS);
	private final SpinnerModel numberOfQuantitativeIons=new SpinnerNumberModel(5, 1, 100, 1);
	private final SpinnerModel minNumOfQuantitativeIons=new SpinnerNumberModel(3, 0, 100, 1);
	private final JComboBox<String> numberOfExtraDecoyLibraries=new JComboBox<String>(NUMBER_OF_EXTRA_DECOY_ITEMS);

	public CASiLParametersPanel() {
		super(new BorderLayout());

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
		options.add(new LabeledComponent("Modification Type", modificationType));
		options.add(new LabeledComponent("Localization Strategy", searchBreadthType));
		options.add(new LabeledComponent("Target/Decoy Approach", numberOfExtraDecoyLibraries));
		options.add(new LabeledComponent("Data Acquisition Type", acquisition));
		options.add(new LabeledComponent("Enzyme", enzyme));
		options.add(new LabeledComponent("Fragmentation", fragType));
		options.add(new LabeledComponent("Precursor Mass Tolerance", precursorTolerance));
		options.add(new LabeledComponent("Fragment Mass Tolerance", fragmentTolerance));
		options.add(new LabeledComponent("Library Mass Tolerance", libraryTolerance));
		options.add(new LabeledComponent("Percolator Version", percolatorVersion));
		options.add(new LabeledComponent("Number of Quantitative Ions", new JSpinner(numberOfQuantitativeIons)));
		options.add(new LabeledComponent("Minimum Number of Quantitative Ions", new JSpinner(minNumOfQuantitativeIons)));
		options.add(new LabeledComponent("Number of Cores", new JSpinner(numberOfJobs)));

		this.add(options, BorderLayout.CENTER);
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
	
	public String getCopy() {
		return copy;
	}
	
	@Override
	public void askForSetupFile() {
		libraryFileChooser.askForFiles();
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.gui.pecan.ParametersPanelInterface#canLoadData()
	 */
	@Override
	public Optional<String> canLoadData() {
		if (libraryFileChooser.getFile()==null) {
			return Optional.of("Please load a library file first!");	
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
		if (libraryFile==null) return null;
		SearchJob job=getJob(diaFile, libraryFile, model, parameters);

		if (job!=null) {
			model.addJob(job);
		}
		return job;
	}

	private static HashMap<File, LibraryInterface> libraries=new HashMap<File, LibraryInterface>();
	static SearchJob getJob(File diaFile, File libraryFile, JobProcessor processor, SearchParameters parameters) {
		File outputFile=new File(diaFile.getAbsolutePath()+CASiLJobData.OUTPUT_FILE_SUFFIX);
		
		LibraryInterface library=libraries.get(libraryFile);
		if (library==null) {
			library=BlibToLibraryConverter.getFile(libraryFile);
			libraries.put(libraryFile, library);
		}
		
		LibraryScoringFactory factory=new EncyclopediaOneScoringFactory(parameters);
		CASiLJobData job=new CASiLJobData(diaFile, library, outputFile, factory);
		return new CASiLJob(processor, job);
	}

	public CASiLSearchParameters getParameters() {
		DataAcquisitionType dataAcquisitionType=DataAcquisitionType.getAcquisitionType((String)acquisition.getSelectedItem());
		DigestionEnzyme digestionEnzyme=DigestionEnzyme.getEnzyme((String)enzyme.getSelectedItem());
		AminoAcidConstants aaConstants=new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());
		FragmentationType fragmentation=FragmentationType.getFragmentationType((String)fragType.getSelectedItem());
		MassTolerance precursorValue=(MassTolerance)precursorTolerance.getSelectedItem();
		MassTolerance fragmentValue=(MassTolerance)fragmentTolerance.getSelectedItem();
		MassTolerance libraryFragmentValue=(MassTolerance)libraryTolerance.getSelectedItem();
		int numberOfJobsValue=((Integer)numberOfJobs.getValue());
		Number value=(Number)precursorWindowWidth.getValue();
		float precursorWindowWidthValue=value==null?-1.0f:value.floatValue();
		boolean isPercolatorTwo=PercolatorExecutor.V2_10.equals(percolatorVersion.getSelectedItem());
		float numberOfExtraDecoyLibrariesValue=NUMBER_OF_EXTRA_DECOY_VALUES[((Integer)numberOfExtraDecoyLibraries.getSelectedIndex())];
		float targetWindowCenter=-1f;
		int numberOfQuantitativeIonsValue=((Integer)numberOfQuantitativeIons.getValue());
		int minNumOfQuantitativeIonsValue=((Integer)minNumOfQuantitativeIons.getValue());
		ScoringBreadthType CASiLSearchBreadthType=(ScoringBreadthType)searchBreadthType.getSelectedItem();
		PeptideModification modification=(PeptideModification)modificationType.getSelectedItem();
		float percolatorThreshold=0.01f;
		CASiLSearchParameters parameters=new CASiLSearchParameters(aaConstants, fragmentation, precursorValue, 0.0, 0.0, fragmentValue, 0.0, libraryFragmentValue, digestionEnzyme, percolatorThreshold, isPercolatorTwo?2:3, dataAcquisitionType, numberOfJobsValue, 25f, targetWindowCenter, precursorWindowWidthValue, numberOfQuantitativeIonsValue, minNumOfQuantitativeIonsValue, modification, CASiLSearchBreadthType, numberOfExtraDecoyLibrariesValue);
		return parameters;
	}
	
	public void setParameters(CASiLSearchParameters params, String libraryFileName) {
		if (libraryFileName!=null) {
			File libraryFile=new File(libraryFileName);
			if (libraryFile.exists()) libraryFileChooser.update(libraryFile);
		}
		
		acquisition.setSelectedItem(DataAcquisitionType.toName(params.getDataAcquisitionType()));
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
		percolatorVersion.setSelectedIndex(params.getPercolatorVersionNumber()==2?1:0);
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
			getParameters().savePreferences(libraryFileChooser.getFile());
		} catch (Exception e) {
			Logger.errorLine("Error writing parameters to disk!");
			Logger.errorException(e);
		}
	}
}
