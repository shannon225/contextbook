package edu.washington.gs.maccoss.encyclopedia;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.UIManager;

import edu.washington.gs.maccoss.encyclopedia.gui.dia.ResultsBrowserPanel;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;

public class DIABrowser {
	public static ImageIcon image=new ImageIcon(DIABrowser.class.getClassLoader().getResource("images/orbi_icon.png"));

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
				System.setProperty("com.apple.mrj.application.apple.menu.about.name", "DIA Browser");
				System.setProperty("apple.laf.useScreenMenuBar", "true");
				break;

			default:
				break;
		}

		EventQueue.invokeLater(new Runnable() {
			@Override
			public void run() {
				final JFrame f=new JFrame("DIA Browser");
				f.setIconImage(image.getImage());
				f.addWindowListener(new WindowAdapter() {
					public void windowClosing(WindowEvent e) {
						System.exit(0);
					}
				});

				f.getContentPane().add(new ResultsBrowserPanel(), BorderLayout.CENTER);

				f.pack();
				f.setSize(new Dimension(2000, 1200));
				f.setVisible(true);
			}
		});

		Logger.logLine("Launching DIA Browser");
	}

}