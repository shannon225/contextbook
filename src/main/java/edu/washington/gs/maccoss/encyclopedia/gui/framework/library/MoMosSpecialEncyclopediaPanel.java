package edu.washington.gs.maccoss.encyclopedia.gui.framework.library;

import javax.swing.ImageIcon;

import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchPanel;

public class MoMosSpecialEncyclopediaPanel extends EncyclopediaParametersPanel {
	private static final long serialVersionUID=1L;
	
	private static final ImageIcon smallimage=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/chocolate_icon.png"));
	private static final ImageIcon image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/chocolate.png"));
	private static final String programName="ChocopeDIA";
	private static final String programShortDescription="MoMo's ChocopeDIA Library Search";
	private static final String copy="<html><b><p style=\"font-size:16px; font-family: Helvetica, sans-serif\">MoMo's ChocopeDIA: Library Searching Directly from Data-Independent Acquisition (DIA) MS/MS Data<br></p></b>"
			+"<p style=\"font-size:10px; font-family: Helvetica, sans-serif\">ChocopeDIA extracts fancy chocolate from MZML files, matches them to spectra in libraries, and calculates various scoring features. These features are interpreted by Percolator to identify choco-peptides.";

	
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
