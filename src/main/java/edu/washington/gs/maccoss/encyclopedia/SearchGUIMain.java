package edu.washington.gs.maccoss.encyclopedia;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.UIManager;

import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchPanel;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;

public class SearchGUIMain {
	public static void main(String[] args) {
		runGUI(false);
	}

	public static void runGUI(boolean pecanpie) {
		try {
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		} catch (Exception e) {
			Logger.errorLine("Error setting look and feel!");
			Logger.errorException(e);
		}

		String shortName;
		String name;
		ImageIcon image;
		if (pecanpie) {
			shortName="PecanPie";
			name="Pecan: Peptide Centric Analysis";
			image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/pecan_icon.png"));
		} else {
			shortName="EncyclopeDIA";
			name="EncyclopeDIA: Peptide Searching for DIA";
			image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/encyclopedia_icon.png"));
		}

		OS os=OSDetector.getOS();
		switch (os) {
			case MAC:
				System.setProperty("com.apple.mrj.application.apple.menu.about.name", shortName);
				System.setProperty("apple.laf.useScreenMenuBar", "true");
				break;

			default:
				break;
		}

		final JFrame f=new JFrame(name);
		f.setIconImage(image.getImage());

		f.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});

		f.getContentPane().add(new SearchPanel(pecanpie), BorderLayout.CENTER);

		f.pack();
		f.setSize(new Dimension(1100, 800));
		f.setVisible(true);

		Logger.logLine(shortName+" Graphical Interface");
	}
}