package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.Color;

import javax.swing.JColorChooser;

public class GUIParameters {
	private static Color BASE_COLOR=new Color(0, 0, 200);

	
	public static Color getBaseColor() {
		return BASE_COLOR;
	}
	public static Color getBrighterColor() {
		int r=Math.min(255, BASE_COLOR.getRed()+100);
		int g=Math.min(255, BASE_COLOR.getGreen()+100);
		int b=Math.min(255, BASE_COLOR.getBlue()+100);
		return new Color(r, g, b, BASE_COLOR.getAlpha());
	}
	public static Color getDarkerColor() {
		int r=Math.max(0, BASE_COLOR.getRed()-100);
		int g=Math.max(0, BASE_COLOR.getGreen()-100);
		int b=Math.max(0, BASE_COLOR.getBlue()-100);
		return new Color(r, g, b, BASE_COLOR.getAlpha());
	}
	
	public static void requestUpdatedColor() {
		BASE_COLOR = JColorChooser.showDialog(null, "Choose a base color", BASE_COLOR);	
	}
}
