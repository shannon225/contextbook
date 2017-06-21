package edu.washington.gs.maccoss.encyclopedia;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.UIManager;

import edu.washington.gs.maccoss.encyclopedia.gui.framework.ParametersPanelInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.MemoryMonitor;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;
import edu.washington.gs.maccoss.encyclopedia.utils.io.Networking;

public class SearchGUIMain {
	public static void main(String[] args) {
		runGUI(ProgramType.Global);
	}

	public static void runGUI(ProgramType program) {
		try {
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		} catch (Exception e) {
			Logger.errorLine("Error setting look and feel!");
			Logger.errorException(e);
		}

		int isOffending=Networking.isOffendingAddress();
		String shortName;
		String name;
		ImageIcon image;
		if (ProgramType.PecanPie==program) {
			shortName="PecanPie";
			name="Pecan: Peptide Centric Analysis";
			image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/pecan_icon.png"));
		} else if (ProgramType.CASiL==program) {
			shortName="CASiL";
			name="CASiL: Chromatogram Aligned Site Localization";
			image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/CASiL_icon.png"));
		} else if (ProgramType.XCorDIA==program) {
			shortName="XCorDIA";
			name="XCorDIA: Peptide Searching with Cross Correlation";
			image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/mike_rotate_icon.png"));
		} else {
			if (isOffending==1) {
				shortName="PoopeDIA";
				name="Lindsay's PoopeDIA: Peptide Searching for DIA";
				image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/mazie_icon.png"));
			} else if (isOffending==2) {
				shortName="ChocopeDIA";
				name="MoMo's ChocopeDIA: Peptide Searching for DIA";
				image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/chocolate2.png"));
			} else {
				shortName="EncyclopeDIA";
				name="EncyclopeDIA: Peptide Searching for DIA";
				image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/encyclopedia_icon.png"));
			}
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

		final SearchPanel panel=new SearchPanel(program);
		f.getContentPane().add(panel, BorderLayout.CENTER);
		f.setJMenuBar(panel.createMenus());

		f.pack();
		f.setSize(new Dimension(1250, 750));
		f.setVisible(true);

		Runtime instance=Runtime.getRuntime();
		long mbOfMemory=(instance.maxMemory()/MemoryMonitor.mb);
		if (mbOfMemory<1000) {
			JOptionPane.showMessageDialog(f, "Warning, you only have "+mbOfMemory+" MB of memory allocated.\nPlease make sure you are running 64-bit Java!", "Warning, Low Memory!", JOptionPane.WARNING_MESSAGE, image);
		}

		Logger.logLine(shortName+" Graphical Interface");
		
		Runtime.getRuntime().addShutdownHook(new Thread() {
			@Override
			public void run() {
				for (ParametersPanelInterface p : panel.getAllTabs()) {
					try {
						p.savePreferences();
					} catch (Exception e) {
						Logger.errorLine("Error writing parameters to disk!");
						Logger.errorException(e);
					}
				}
			}
		});
	}
}