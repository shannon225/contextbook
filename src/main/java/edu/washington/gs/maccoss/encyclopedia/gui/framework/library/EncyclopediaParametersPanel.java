package edu.washington.gs.maccoss.encyclopedia.gui.framework.library;

import java.awt.BorderLayout;
import java.awt.Color;
import java.io.File;
import java.text.NumberFormat;
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

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
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

public class EncyclopediaParametersPanel extends JPanel implements ParametersPanelInterface {
	
	private static final long serialVersionUID=1L;
	private static final int numberOfCores=Runtime.getRuntime().availableProcessors();
	private static final String PHOSPHOPROTEOME="Phosphoproteome";
	private static final String[] NUMBER_OF_EXTRA_DECOY_ITEMS=new String[] {"Normal Target/Decoy", "+10% Extra Decoys", "+20% Extra Decoys", "+50% Extra Decoys", "+100% Extra Decoys (2x Time)"};
	private static final float[] NUMBER_OF_EXTRA_DECOY_VALUES=new float[] {0.0f, 0.1f, 0.2f, 0.5f, 1.0f};
	
	private static final ImageIcon smallimage=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/encyclopedia_small_icon.png"));
	private static final ImageIcon image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/encyclopedia_icon.png"));
	private static final String programName="EncyclopeDIA";
	private static final String programShortDescription="EncyclopeDIA Library Search";
	private static final String copy="<html><b><p style=\"font-size:16px; font-family: Helvetica, sans-serif\">EncyclopeDIA: Library Searching Directly from Data-Independent Acquisition (DIA) MS/MS Data<br></p></b>"
			+ "<p style=\"font-size:10px; font-family: Helvetica, sans-serif\">EncyclopeDIA extracts peptide fragmentation chromatograms from MZML files, matches them to spectra in libraries, and calculates various scoring features. These features are interpreted by Percolator to identify peptides.";
	
	private final FileChooserPanel libraryFileChooser;
	private final JComboBox<String> acquisition=new JComboBox<String>(new String[] {DataAcquisitionType.toName(DataAcquisitionType.OVERLAPPING_DIA), DataAcquisitionType.toName(DataAcquisitionType.DIA), DataAcquisitionType.toName(DataAcquisitionType.DDA)});
	private final JComboBox<String> enzyme=new JComboBox<String>(new String[] {"Trypsin", "Lys-C", "Lys-N", "Arg-C", "CNBr", "Chymotrypsin", "PepsinA"});
	private final JComboBox<String> fragType=new JComboBox<String>(new String[] {FragmentationType.toName(FragmentationType.CID), FragmentationType.toName(FragmentationType.YONLY), FragmentationType.toName(FragmentationType.ETD)});
	private final JComboBox<String> proteomeType=new JComboBox<String>(new String[] {"Standard Proteome", PHOSPHOPROTEOME});

	private final JFormattedTextField precursorWindowWidth=new JFormattedTextField(NumberFormat.getNumberInstance()); // not displayed anymore

	private final SpinnerModel precursorPPM=new SpinnerNumberModel(10, 1, 1000, 1);
	private final SpinnerModel fragmentPPM=new SpinnerNumberModel(10, 1, 1000, 1);
	private final SpinnerModel libraryFragmentPPM=new SpinnerNumberModel(10, 1, 1000, 1);
	private final SpinnerModel numberOfJobs=new SpinnerNumberModel(numberOfCores, 1, numberOfCores, 1);
	private final SpinnerModel numberOfQuantitativeIons=new SpinnerNumberModel(5, 1, 100, 1);
	private final JComboBox<String> numberOfExtraDecoyLibraries=new JComboBox<String>(NUMBER_OF_EXTRA_DECOY_ITEMS);

	public EncyclopediaParametersPanel() {
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
		
		libraryFileChooser=new FileChooserPanel(null, "Library", new SimpleFilenameFilter(".elib"), true);
		options.add(libraryFileChooser);
		options.add(new LabeledComponent("Target/Decoy Approach", numberOfExtraDecoyLibraries));
		options.add(new LabeledComponent("Data Acquisition Type", acquisition));
		options.add(new LabeledComponent("Enzyme", enzyme));
		options.add(new LabeledComponent("Fragmentation", fragType));
		options.add(new LabeledComponent("Proteome Type", proteomeType));
		options.add(new LabeledComponent("Precursor (PPM)", new JSpinner(precursorPPM)));
		options.add(new LabeledComponent("Fragment (PPM)", new JSpinner(fragmentPPM)));
		options.add(new LabeledComponent("Library Fragment (PPM)", new JSpinner(libraryFragmentPPM)));
		options.add(new LabeledComponent("Number of Quantitative Ions", new JSpinner(numberOfQuantitativeIons)));
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

	static SearchJob getJob(File diaFile, File libraryFile, JobProcessor processor, SearchParameters parameters) {
		File outputFile=new File(diaFile.getAbsolutePath()+EncyclopediaJobData.OUTPUT_FILE_SUFFIX);
		
		LibraryInterface library=BlibToLibraryConverter.getFile(libraryFile);
		
		LibraryScoringFactory factory=new EncyclopediaOneScoringFactory(parameters);
		EncyclopediaJobData job=new EncyclopediaJobData(diaFile, library, outputFile, factory);
		return new EncyclopediaJob(processor, job);
	}

	public SearchParameters getParameters() {
		DataAcquisitionType dataAcquisitionType=DataAcquisitionType.getAcquisitionType((String)acquisition.getSelectedItem());
		DigestionEnzyme digestionEnzyme=DigestionEnzyme.getEnzyme((String)enzyme.getSelectedItem());
		AminoAcidConstants aaConstants=new AminoAcidConstants(new TCharDoubleHashMap());
		FragmentationType fragmentation=FragmentationType.getFragmentationType((String)fragType.getSelectedItem());
		float precursorPPMValue=((Number)precursorPPM.getValue()).floatValue();
		float fragmentPPMValue=((Number)fragmentPPM.getValue()).floatValue();
		float libraryFragmentPPMValue=((Number)libraryFragmentPPM.getValue()).floatValue();
		int numberOfJobsValue=((Integer)numberOfJobs.getValue());
		Number value=(Number)precursorWindowWidth.getValue();
		float precursorWindowWidthValue=value==null?-1.0f:value.floatValue();
		boolean isPhospho=PHOSPHOPROTEOME.equals(proteomeType.getSelectedItem());
		float numberOfExtraDecoyLibrariesValue=NUMBER_OF_EXTRA_DECOY_VALUES[((Integer)numberOfExtraDecoyLibraries.getSelectedIndex())];
		float targetWindowCenter=-1f;
		int numberOfQuantitativeIonsValue=((Integer)numberOfQuantitativeIons.getValue());
		SearchParameters parameters=new SearchParameters(aaConstants, fragmentation, new MassTolerance(precursorPPMValue, MassErrorUnitType.PPM), 0.0, new MassTolerance(fragmentPPMValue, MassErrorUnitType.PPM), 0.0, new MassTolerance(libraryFragmentPPMValue, MassErrorUnitType.PPM), digestionEnzyme, 0.01f, null, dataAcquisitionType, numberOfJobsValue, 25f, targetWindowCenter, precursorWindowWidthValue, numberOfQuantitativeIonsValue, isPhospho, numberOfExtraDecoyLibrariesValue);
		return parameters;
	}
	
	public void setParameters(SearchParameters params, String libraryFileName) {
		if (libraryFileName!=null) {
			File libraryFile=new File(libraryFileName);
			if (libraryFile.exists()) libraryFileChooser.update(libraryFile);
		}
		
		acquisition.setSelectedItem(DataAcquisitionType.toName(params.getDataAcquisitionType()));
		enzyme.setSelectedItem(params.getEnzyme().getName());
		fragType.setSelectedItem(FragmentationType.toName(params.getFragType()));
		precursorPPM.setValue((int)params.getPrecursorTolerance().getPpmTolerance());
		fragmentPPM.setValue((int)params.getFragmentTolerance().getPpmTolerance());
		libraryFragmentPPM.setValue((int)params.getLibraryFragmentTolerance().getPpmTolerance());
		numberOfJobs.setValue(params.getNumberOfThreadsUsed());
		if (params.getPrecursorWindowSize()>0) {
			precursorWindowWidth.setValue(params.getPrecursorWindowSize());
		} else {
			precursorWindowWidth.setValue(-1);
		}
		proteomeType.setSelectedIndex(params.isRunPhosphoLocalization()?1:0);
		int index=Arrays.binarySearch(NUMBER_OF_EXTRA_DECOY_VALUES, params.getNumberOfExtraDecoyLibrariesSearched());
		if (index>=0) {
			numberOfExtraDecoyLibraries.setSelectedIndex(index);
		}
		numberOfQuantitativeIons.setValue(params.getNumberOfQuantitativePeaks());
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
