package edu.washington.gs.maccoss.encyclopedia.gui.framework.library;

import javax.swing.ImageIcon;

import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchPanel;

public class AustinsSpecialEncyclopediaPanel extends EncyclopediaParametersPanel {
	private static final long serialVersionUID=1L;
	
	private static final ImageIcon smallimage=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/austin_small_icon.png"));
	private static final ImageIcon image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/austin_icon.png"));
	private static final String programName="StupiDIA";
	private static final String programShortDescription="Austin's Stupid Library Search";
	private static final String copy="<html><b><p style=\"font-size:16px; font-family: Helvetica, sans-serif\">Austin's StupiDIA Library Search<br></p></b>"
			+"<p style=\"font-size:10px; font-family: Helvetica, sans-serif\">StupiDIA extracts peptide fragmentation chromatograms from MZML files, matches them to spectra in libraries, and calculates various scoring features. These features are interpreted by Percolator to identify peptides.";

	public AustinsSpecialEncyclopediaPanel(SearchPanel searchPanel) {
		super(searchPanel);
	}

	public String getProgramName() {
		return programName;
	}
	
	@Override
	public String getProgramShortDescription() {
		return programShortDescription;
	}
	
	@Override
	public ImageIcon getSmallImage() {
		return smallimage;
	}

	@Override
	public ImageIcon getImage() {
		return image;
	}

	@Override
	public String getCopy() {
		return copy;
	}
}
