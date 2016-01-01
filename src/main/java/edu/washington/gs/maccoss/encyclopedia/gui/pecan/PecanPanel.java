package edu.washington.gs.maccoss.encyclopedia.gui.pecan;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FileDialog;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.table.TableColumn;

import com.google.common.base.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LabeledComponent;
import edu.washington.gs.maccoss.encyclopedia.gui.general.ProgressRenderer;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;

public class PecanPanel extends JPanel {
	private static final long serialVersionUID=1L;

	public static final String copy="<html><b><p style=\"font-size:20px; font-family: Helvetica, sans-serif\">PECAN: Peptide Detection Directly from Data-Independent Acquisition (DIA) MS/MS Data<br><br></p>"
			+ "<p style=\"font-size:12px; font-family: Helvetica, sans-serif\">PECAN extracts peptide fragmentation chromatograms from MZML files, assigns peaks, and calculates various peak features. These features are interpreted by Percolator to identify peptides.";
	
	private final FileChooserPanel backgroundFasta;
	private final JComboBox<String> enzyme=new JComboBox<String>(new String[] {"Trypsin", "Lys-C", "Lys-N", "Arg-C", "CNBr", "Chymotrypsin", "PepsinA"});
	private final JComboBox<String> fixed=new JComboBox<String>(new String[] {"C+57 (Carbamidomethyl)", "C+58 (Carboxymethyl)", "C+46 (MMTS)", "None"});
	private final JComboBox<String> fragType=new JComboBox<String>(new String[] {"HCD (Y-Only)", "CID (B/Y)", "ETD (C/Z/Z+1)"});

	private final SpinnerModel precursorPPM = new SpinnerNumberModel(10, 1, 1000, 1);
	private final SpinnerModel fragmentPPM = new SpinnerNumberModel(10, 1, 1000, 1);
	private final SpinnerModel minCharge = new SpinnerNumberModel(2, 1, 2, 1);
	private final SpinnerModel maxCharge = new SpinnerNumberModel(3, 2, 4, 1);
	private final SpinnerModel maxMissedCleavage = new SpinnerNumberModel(1, 0, 3, 1);
	
	PecanFileProcessorModel pecanModel=new PecanFileProcessorModel();
	
	public PecanPanel() {
		super(new BorderLayout());
	    
		setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
		setBackground(Color.white);

		JPanel top=new JPanel(new BorderLayout());
		ImageIcon image=new ImageIcon(this.getClass().getClassLoader().getResource("images/pecan.png"));
		top.add(new JLabel(image), BorderLayout.WEST);
		JEditorPane editor=new JEditorPane("text/html", copy);
		editor.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
		top.add(editor, BorderLayout.CENTER);
		top.setBackground(Color.white);

		this.add(top, BorderLayout.NORTH);

		JSplitPane split=new JSplitPane();

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(new LabeledComponent("<p style=\"font-size:12px; font-family: Helvetica, sans-serif\"><b>Parameters", new JLabel()));
		
		backgroundFasta=new FileChooserPanel(null, "Background", new SimpleFilenameFilter(".fas", ".fasta"));
		options.add(backgroundFasta);
		options.add(new LabeledComponent("Enzyme", enzyme));
		options.add(new LabeledComponent("Fixed", fixed));
		options.add(new LabeledComponent("Fragmentation", fragType));
		options.add(new LabeledComponent("Precursor (PPM)", new JSpinner(precursorPPM)));
		options.add(new LabeledComponent("Fragment (PPM)", new JSpinner(fragmentPPM)));
		options.add(new LabeledComponent("Maximum Missed Cleavage", new JSpinner(maxMissedCleavage)));

		JPanel chargeRange=new JPanel(new FlowLayout());
		chargeRange.setOpaque(true);
		chargeRange.setBackground(Color.white);
		chargeRange.add(new JSpinner(minCharge));
		chargeRange.add(new JLabel("<html><p style=\"font-size:10px; font-family: Helvetica, sans-serif\"> to "));
		chargeRange.add(new JSpinner(maxCharge));
		options.add(new LabeledComponent("Charge range", chargeRange));

		
		JPanel optionsWrapper=new JPanel(new BorderLayout());
		optionsWrapper.setOpaque(true);
		optionsWrapper.setBackground(Color.white);
		optionsWrapper.add(options, BorderLayout.NORTH);
		
		split.setLeftComponent(optionsWrapper);

		JPanel files=new JPanel(new BorderLayout());
		JButton chooseFile=new JButton("Add");
		chooseFile.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				JFrame frame = (JFrame)SwingUtilities.getRoot(PecanPanel.this);
				FileDialog dialog=new FileDialog(frame, "Select a MZML file", FileDialog.LOAD);
				dialog.setFilenameFilter(new SimpleFilenameFilter(".mzml", ".dia"));
				dialog.setVisible(true);
				if (dialog.getFiles()!=null) {
					for (File file : dialog.getFiles()) {
						PecanJob job=getJob(file);
						if (job!=null) {
							pecanModel.addJob(job);
						}
					}
				}
			}
		});
		
		files.add(new JLabel("<html><p style=\"font-size:12px; font-family: Helvetica, sans-serif\"><b>MZML Files: "), BorderLayout.WEST);
		files.add(chooseFile, BorderLayout.EAST);
		
		JPanel filesWrapper=new JPanel(new BorderLayout());
		filesWrapper.setOpaque(true);
		filesWrapper.setBackground(Color.white);
		filesWrapper.add(files, BorderLayout.NORTH);
		
		JTable table=new JTable(pecanModel);
        TableColumn column = table.getColumnModel().getColumn(1);
        column.setCellRenderer(new ProgressRenderer());
        
		filesWrapper.add(table, BorderLayout.CENTER);
		
		split.setRightComponent(filesWrapper);
		split.setDividerLocation(300);

		this.add(split, BorderLayout.CENTER);
	}
	
	public PecanJob getJob(File diaFile) {
		File fastaFile=backgroundFasta.getFile();
		if (fastaFile==null) return null;
		
		File outputFile=new File(diaFile.getAbsolutePath()+".pecan.txt");
		File featureFile=new File(outputFile.getAbsolutePath()+".features.txt");
		
		ArrayList<FastaEntry> targets=null;
		DigestionEnzyme digestionEnzyme=DigestionEnzyme.getEnzyme((String)enzyme.getSelectedItem());
		AminoAcidConstants aaConstants=AminoAcidConstants.getConstants((String)fixed.getSelectedItem());
		FragmentationType fragmentation=FragmentationType.getFragmentationType((String)fragType.getSelectedItem());
		float precursorPPMValue=((Integer)precursorPPM.getValue()).floatValue();
		float fragmentPPMValue=((Integer)fragmentPPM.getValue()).floatValue();
		byte minChargeValue=((Integer)minCharge.getValue()).byteValue();
		byte maxChargeValue=((Integer)maxCharge.getValue()).byteValue();
		byte maxMissedCleavageValue=((Integer)maxMissedCleavage.getValue()).byteValue();
		SearchParameters parameters=new SearchParameters(aaConstants, fragmentation, new MassTolerance(precursorPPMValue), new MassTolerance(fragmentPPMValue), digestionEnzyme,
				maxMissedCleavageValue, minChargeValue, maxChargeValue);
		
		PecanScoringFactory factory=new PecanOneScoringFactory(parameters, featureFile);
		
		return new PecanJob(pecanModel, Optional.fromNullable(targets), diaFile, fastaFile, featureFile, outputFile, factory);
	}
}
