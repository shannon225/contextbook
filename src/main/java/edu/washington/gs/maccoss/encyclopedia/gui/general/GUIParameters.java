package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.Color;
import java.io.IOException;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

import javax.swing.JColorChooser;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;

public class GUIParameters {
	private static Nothing initialized=null;
	private static Color BASE_COLOR=new Color(0, 0, 200);

	private synchronized static void checkInit() {
		if (initialized==null) {
			BASE_COLOR=new Color(0, 0, 200);
		
			try {
				readPreferences();
			} catch (Exception e) {
				Logger.errorLine("Error reading GUI preferences!");
			} finally {
				initialized=Nothing.NOTHING;
			}
		}
	}
	
	public static void savePreferences() throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("encyclopedia_gui");
		prefs.put("BASE_COLOR", Integer.toString(BASE_COLOR.getRGB()));
		prefs.flush();
	}
	
	public static void readPreferences() throws IOException,BackingStoreException {
		Preferences prefs=Preferences.userRoot().node("encyclopedia_gui");
		String value=prefs.get("BASE_COLOR", Integer.toString(BASE_COLOR.getRGB()));
		BASE_COLOR=new Color(Integer.parseInt(value));
	}
	
	public static Color getBaseColor() {
		checkInit();
		
		return BASE_COLOR;
	}
	
	public static Color getBaseColor(int alpha) {
		checkInit();

		return new Color(BASE_COLOR.getRed(), BASE_COLOR.getGreen(), BASE_COLOR.getBlue(), alpha);
	}
	
	public static Color getBrighterColor() {
		return getBrighterColor(BASE_COLOR.getAlpha());
	}

	public static Color getBrighterColor(int alpha) {
		checkInit();
		
		int r=Math.min(255, BASE_COLOR.getRed()+100);
		int g=Math.min(255, BASE_COLOR.getGreen()+100);
		int b=Math.min(255, BASE_COLOR.getBlue()+100);
		return new Color(r, g, b, alpha);
	}
	
	public static Color getDarkerColor() {
		return getDarkerColor(BASE_COLOR.getAlpha());
	}
	public static Color getDarkerColor(int alpha) {
		checkInit();
		
		int r=Math.max(0, BASE_COLOR.getRed()-100);
		int g=Math.max(0, BASE_COLOR.getGreen()-100);
		int b=Math.max(0, BASE_COLOR.getBlue()-100);
		return new Color(r, g, b, alpha);
	}
	
	public static void requestUpdatedColor() {
		checkInit();
		
		Color c=JColorChooser.showDialog(null, "Choose a base color", BASE_COLOR);
		if (c!=null) {
			BASE_COLOR = c;
			try {
				savePreferences();
			} catch (Exception e) {
				Logger.errorLine("Error saving GUI preferences!");
			}
		}
	}
}
