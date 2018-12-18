package edu.washington.gs.maccoss.encyclopedia.utils;

import java.awt.Desktop;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.net.URL;

import javax.swing.ImageIcon;
import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;

import com.brsanthu.googleanalytics.GoogleAnalytics;

import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.utils.io.Version;

public class VersioningDetector {
	private static final String ENCYCLOPEDIA_URL= "https://bitbucket.org/searleb/encyclopedia/";

	public static boolean checkVersionCLI(ProgramType program) {
		return checkVersionGUI(program, null);
	}
	public static boolean checkVersionGUI(ProgramType program, JFrame frame) {
		try {
			Version localVersion=program.getVersion();
			if (localVersion.amIAbove(new Version(0, 0, 0))) { // otherwise is a development version
				GoogleAnalytics ga=GoogleAnalytics.builder().withTrackingId("UA-131121966-1").withAppName(program.toString()).withAppVersion(program.toString()+" "+localVersion.toString()).build();
				ga.screenView(program.toString(), localVersion.toString()).send();
				URL url = new URL(ENCYCLOPEDIA_URL+"downloads/current_version.txt");
				BufferedReader in = new BufferedReader(new InputStreamReader(url.openStream()));
				String downloadedVersion=in.readLine();
				Version onlineVersion=new Version(downloadedVersion);
				if (onlineVersion.amIAbove(localVersion)) {
					Logger.errorLine("I found a new version on the web! You're running "+localVersion+" and you should update to "+onlineVersion);
					if (frame!=null) {
						announceUpdateNotice(program, frame, localVersion, onlineVersion);
					}
					return false; // should update
				}
			}
		} catch (Exception e) {
			Logger.logLine("Sorry, I ran into an error checking for new versions. You should look online at ["+ENCYCLOPEDIA_URL+"] if you need to update your version!");
		}
		return true;
	}

	private static void announceUpdateNotice(ProgramType program, JFrame frame, Version localVersion, Version onlineVersion) {
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
	    
		JOptionPane.showMessageDialog(frame, message, "New version of " + program.toString(),
				JOptionPane.INFORMATION_MESSAGE, icon);
	}
}
