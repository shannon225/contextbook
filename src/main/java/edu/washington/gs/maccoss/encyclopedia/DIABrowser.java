package edu.washington.gs.maccoss.encyclopedia;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.HashMap;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.UIManager;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.dia.ResultsBrowserPanel;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;
import edu.washington.gs.maccoss.encyclopedia.utils.VersioningDetector;

public class DIABrowser {
	public static ImageIcon image=new ImageIcon(DIABrowser.class.getClassLoader().getResource("images/orbi_icon.png"));

	public static void main(String[] args) {
		VersioningDetector.checkVersionCLI(ProgramType.EncyclopeDIA);
		
		//final File dia=new File("/Users/searleb/Documents/school/localization_manuscript/elibs/110515_bcs_hela_phospho_starved_20mz_500_900.dia");
		//final File library=new File("/Users/searleb/Documents/school/localization_manuscript/elibs/110515_bcs_hela_phospho_starved_20mz_500_900.dia.elib");
		//final File dia=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/dec2015_phospho/110515_bcs_hela_phospho_starved_20mz_500_900.dia");
		//final File library=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/dec2015_phospho/110515_bcs_hela_phospho_starved_20mz_500_900.dia.elib");

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
				
				HashMap<String, String> map=SearchParameterParser.getDefaultParameters();
				map.put("-deconvoluteOverlappingWindows", "true");
				map.put("-fixed", "");
				SearchParameters parameters=SearchParameterParser.parseParameters(map);

				ResultsBrowserPanel panel=new ResultsBrowserPanel(parameters);
				//panel.updateTable(library);
				//panel.updateRaw(dia);
				f.getContentPane().add(panel, BorderLayout.CENTER);

				f.pack();
				f.setSize(new Dimension(1900, 1030)); // for 1920x1080
				f.setVisible(true);
			}
		});

		Logger.logLine("Launching DIA Browser");
	}

}