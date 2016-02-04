package edu.washington.gs.maccoss.encyclopedia;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.JFrame;
import javax.swing.UIManager;

import edu.washington.gs.maccoss.encyclopedia.gui.dia.DIABrowserPanel;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;

public class DIABrowser {
	//public static ImageIcon image=new ImageIcon(DIABrowser.class.getClassLoader().getResource("images/orbi.png"));
	
	public static void main(String[] args) {
		try {
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		} catch (Exception e) {
			Logger.errorLine("Error setting look and feel!");
			Logger.errorException(e);
		}
		OS os=OSDetector.getOS();
		switch (os) {
			case MAC:
				System.setProperty("com.apple.mrj.application.apple.menu.about.name", "DIABrowser");
				System.setProperty("apple.laf.useScreenMenuBar", "true");
				break;

			default:
				break;
		}
		
		final JFrame f=new JFrame("DIA Browser");
		//f.setIconImage(image.getImage());
		f.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});

		f.getContentPane().add(new DIABrowserPanel(), BorderLayout.CENTER);
		
		f.pack();
		f.setSize(new Dimension(1100, 800));
		f.setVisible(true);

		Logger.logLine("Launching DIA Browser");
	}

}