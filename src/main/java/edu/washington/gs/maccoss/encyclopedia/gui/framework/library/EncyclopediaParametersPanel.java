package edu.washington.gs.maccoss.encyclopedia.gui.framework.library;

import java.awt.BorderLayout;
import java.awt.Color;
import java.io.File;
import java.text.NumberFormat;
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

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DataAcquisitionType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
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
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import gnu.trove.map.hash.TCharFloatHashMap;

public class EncyclopediaParametersPanel extends JPanel implements ParametersPanelInterface {
	private static final long serialVersionUID=1L;
	private static final int numberOfCores=Runtime.getRuntime().availableProcessors();
	private static final ImageIcon smallimage=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/encyclopedia_small_icon.png"));
	private static final ImageIcon image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/encyclopedia_icon.png"));
	private static final String programName="EncyclopeDIA";
	private static final String programShortDescription="EncyclopeDIA Library Search";
	private static final String copy="<html><b><p style=\"font-size:16px; font-family: Helvetica, sans-serif\">EncyclopeDIA: Library Searching Directly from Data-Independent Acquisition (DIA) MS/MS Data<br></p></b>"
			+ "<p style=\"font-size:10px; font-family: Helvetica, sans-serif\">EncyclopeDIA extracts peptide fragmentation chromatograms from MZML files, matches them to spectra in libraries, and calculates various scoring features. These features are interpreted by Percolator to identify peptides.";
	
	private final FileChooserPanel libraryFileChooser;
	private final JComboBox<String> acquisition=new JComboBox<String>(new String[] {"Overlapping DIA", "Non-Overlapping DIA", "DDA"});
	private final JComboBox<String> enzyme=new JComboBox<String>(new String[] {"Trypsin", "Lys-C", "Lys-N", "Arg-C", "CNBr", "Chymotrypsin", "PepsinA"});
	private final JComboBox<String> fragType=new JComboBox<String>(new String[] {"CID (B/Y)", "HCD (Y-Only)", "ETD (C/Z/Z+1)"});

	private final JFormattedTextField precursorWindowWidth=new JFormattedTextField(NumberFormat.getNumberInstance());

	private final SpinnerModel precursorPPM = new SpinnerNumberModel(10, 1, 1000, 1);
	private final SpinnerModel fragmentPPM = new SpinnerNumberModel(10, 1, 1000, 1);
	private final SpinnerModel numberOfJobs = new SpinnerNumberModel(numberOfCores, 1, numberOfCores, 1);

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
		options.add(new LabeledComponent("Data Acquisition Type", acquisition));
		options.add(new LabeledComponent("Precursor Window Width (blank=extract from file)", precursorWindowWidth));
		options.add(new LabeledComponent("Enzyme", enzyme));
		options.add(new LabeledComponent("Fragmentation", fragType));
		options.add(new LabeledComponent("Precursor (PPM)", new JSpinner(precursorPPM)));
		options.add(new LabeledComponent("Fragment (PPM)", new JSpinner(fragmentPPM)));
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
		File outputFile=new File(diaFile.getAbsolutePath()+".encyclopedia.txt");
		File featureFile=new File(outputFile.getAbsolutePath()+".features.txt");
		
		LibraryInterface library=BlibToLibraryConverter.getFile(libraryFile);
		
		LibraryScoringFactory factory=new EncyclopediaOneScoringFactory(parameters);
		EncyclopediaJobData job=new EncyclopediaJobData(diaFile, library, featureFile, outputFile, factory);
		return new EncyclopediaJob(processor, job);
	}

	private SearchParameters getParameters() {
		DataAcquisitionType dataAcquisitionType=DataAcquisitionType.getAcquisitionType((String)acquisition.getSelectedItem());
		DigestionEnzyme digestionEnzyme=DigestionEnzyme.getEnzyme((String)enzyme.getSelectedItem());
		AminoAcidConstants aaConstants=new AminoAcidConstants(new TCharFloatHashMap());
		FragmentationType fragmentation=FragmentationType.getFragmentationType((String)fragType.getSelectedItem());
		float precursorPPMValue=((Integer)precursorPPM.getValue()).floatValue();
		float fragmentPPMValue=((Integer)fragmentPPM.getValue()).floatValue();
		int numberOfJobsValue=((Integer)numberOfJobs.getValue());
		Number value=(Number)precursorWindowWidth.getValue();
		float precursorWindowWidthValue=value==null?-1.0f:value.floatValue();
		SearchParameters parameters=new SearchParameters(aaConstants, fragmentation, new MassTolerance(precursorPPMValue), new MassTolerance(fragmentPPMValue), digestionEnzyme, 0.01f, null, dataAcquisitionType, numberOfJobsValue, 25f, -1f, precursorWindowWidthValue, true, true);
		return parameters;
	}
}
