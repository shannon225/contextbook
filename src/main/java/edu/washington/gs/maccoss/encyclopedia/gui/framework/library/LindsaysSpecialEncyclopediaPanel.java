package edu.washington.gs.maccoss.encyclopedia.gui.framework.library;

import javax.swing.ImageIcon;

import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchPanel;

public class LindsaysSpecialEncyclopediaPanel extends EncyclopediaParametersPanel {
	private static final long serialVersionUID=1L;
	
	private static final ImageIcon smallimage=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/mazie_small_icon.png"));
	private static final ImageIcon image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/mazie_icon.png"));
	private static final String programName="MaiziepeDIA";
	private static final String programShortDescription="Lindsay's MaiziepeDIA Library Search";
	private static final String copy="<html><b><p style=\"font-size:16px; font-family: Helvetica, sans-serif\">Lindsay's MaiziepeDIA: Now powered by chilly pup!<br></p></b>"
			+"<p style=\"font-size:10px; font-family: Helvetica, sans-serif\">MaiziepeDIA extracts peptide fragmentation chromatograms from MZML files, matches them to spectra in libraries, and calculates various scoring features. These features are interpreted by Percolator to identify peptides.";

	public LindsaysSpecialEncyclopediaPanel(SearchPanel searchPanel) {
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
