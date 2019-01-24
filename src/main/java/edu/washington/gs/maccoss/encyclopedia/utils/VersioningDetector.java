package edu.washington.gs.maccoss.encyclopedia.utils;

import java.awt.Component;
import java.awt.Desktop;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

import javax.swing.ImageIcon;
import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;

import com.brsanthu.googleanalytics.GoogleAnalytics;

import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.utils.io.Networking;
import edu.washington.gs.maccoss.encyclopedia.utils.io.Version;

public class VersioningDetector {
	private static final String ENCYCLOPEDIA_URL= "https://bitbucket.org/searleb/encyclopedia/";
	
	public static void main(String[] args) {
		VersioningDetector.checkVersionGUI(new Version(0,7,1), new Version(0,7,2), ProgramType.EncyclopeDIA, new JFrame(), false, true);
	}

	public static boolean checkVersionCLI(ProgramType program) {
		return checkVersionGUI(program, null);
	}
	public static boolean checkVersionGUI(ProgramType program, Component frame) {
		return checkVersionGUI(program, frame, false, false);
	}
	public static boolean checkVersionGUI(ProgramType program, Component frame, boolean ignoreDev, boolean forceReportToUser) {
		Version localVersion=program.getVersion();
		return checkVersionGUI(localVersion, program, frame, ignoreDev, forceReportToUser);
	}
	public static boolean checkVersionGUI(Version localVersion, ProgramType program, Component frame, boolean ignoreDev, boolean forceReportToUser) {
		if (ignoreDev||localVersion.amIAbove(new Version(0, 0, 0))) { // otherwise is a development version
			try {
				Version onlineVersion=getOnlineVersion(localVersion, program);
				return checkVersionGUI(localVersion, onlineVersion, program, frame, ignoreDev, forceReportToUser);
			} catch (Exception e) {
				Logger.logLine("Sorry, I ran into an error checking for new versions. You should look online at ["+ENCYCLOPEDIA_URL+"] if you need to update your version!");
				return false; // should update
			}
		} else {
			return true;
		}
	}
	public static boolean checkVersionGUI(Version localVersion, Version onlineVersion, ProgramType program, Component frame, boolean ignoreDev, boolean forceReportToUser) {
		if (ignoreDev||localVersion.amIAbove(new Version(0, 0, 0))) { // otherwise is a development version
			if (onlineVersion.amIAbove(localVersion)) {
				Version dontCheckBelow=readStopCheckingUntilPreference();
				Logger.errorLine("I found a new version on the web! You're running "+localVersion+" and you should update to "+onlineVersion);
				if (onlineVersion.amIAbove(dontCheckBelow)&&frame!=null) {
					announceUpdateNotice(program, frame, localVersion, onlineVersion);
				}
				return false; // should update
			}
		}
		return true;
	}

	private static Version getOnlineVersion(Version localVersion, ProgramType program) throws MalformedURLException, IOException {
		GoogleAnalytics ga=GoogleAnalytics.builder().withTrackingId("UA-131121966-1").withAppName(program.toString()).withAppVersion(program.toString()+" "+localVersion.toString()).build();
		ga.screenView(program.toString(), localVersion.toString()).userId(Networking.getUserID()).sendAsync();
		URL url = new URL(ENCYCLOPEDIA_URL+"downloads/current_version.txt");
		BufferedReader in = new BufferedReader(new InputStreamReader(url.openStream()));
		String downloadedVersion=in.readLine();
		Version onlineVersion=new Version(downloadedVersion);
		return onlineVersion;
	}
	

	
	public static void stopCheckingUntilNewVersion(Version v) {
		Preferences prefs=Preferences.userRoot().node("encyclopedia_versioning");
		prefs.put("version", v.toString());
		
		try {
			prefs.flush();
		} catch (BackingStoreException bse) {
			Logger.errorLine("Ran into error saving ignore version to disk, continuing without saving.");
			Logger.errorException(bse);
		}
	}
	
	public static Version readStopCheckingUntilPreference() {
		Preferences prefs=Preferences.userRoot().node("encyclopedia_versioning");
		String v=prefs.get("version", null);
		if (v==null) {
			return new Version(0, 0, 0);
		}
		return new Version(v);
	}

	private static void announceUpdateNotice(ProgramType program, Component frame, Version localVersion, Version onlineVersion) {
		ImageIcon icon=new ImageIcon(VersioningDetector.class.getClassLoader().getResource("images/encyclopedia_icon.png"));

	    JLabel label = new JLabel();
	    
	    JEditorPane message = new JEditorPane("text/html", "<html><p style=\"font-size:12px; font-family: Helvetica, sans-serif\">There's a new version of <b>"+program.toString()+"</b> on the web!<br>The new version is "+onlineVersion+" and you are running "+localVersion+".<br>You can update at <a href="+ENCYCLOPEDIA_URL+">"+ENCYCLOPEDIA_URL+"</a></html>");
	    message.addHyperlinkListener(new HyperlinkListener()
	    {
	        @Override
	        public void hyperlinkUpdate(HyperlinkEvent e)
	        {
	            if (e.getEventType().equals(HyperlinkEvent.EventType.ACTIVATED)) {
		            if (Desktop.isDesktopSupported()) {
		                try {
		                    Desktop.getDesktop().browse(e.getURL().toURI());
		                } catch (IOException e1) {
		                    Logger.errorException(e1);
		                } catch (URISyntaxException e1) {
		                    Logger.errorException(e1);
		                }
		            }
	            }
	        }
	    });
	    message.setEditable(false);
	    message.setBackground(label.getBackground());
		
		Object[] options = { "Ok", "Ignore this version" };
		int ret=JOptionPane.showOptionDialog(frame, message, "New version of " + program.toString(),
		    JOptionPane.DEFAULT_OPTION, JOptionPane.PLAIN_MESSAGE, icon, 
		    options, options[0]);
		if (ret==1) {
			Logger.logLine("User opted to ignore "+onlineVersion);
			stopCheckingUntilNewVersion(onlineVersion);
		}
	}
}
