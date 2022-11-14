package edu.washington.gs.maccoss.encyclopedia;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.HashMap;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.UIManager;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.ParametersPanelInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.MemoryMonitor;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;
import edu.washington.gs.maccoss.encyclopedia.utils.VersioningDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.io.Networking;

public class SearchGUIMain {
	public static void main(String[] args) {
		HashMap<String, String> arguments=CommandLineParser.parseArguments(args);
		runGUI(ProgramType.Global, arguments.containsKey(SearchParameters.ENABLE_ADVANCED_OPTIONS));
	}
	
	public static JFrame runGUI(ProgramType program) {
		return runGUI(program, false);
	}

	public static JFrame runGUI(ProgramType program, boolean enableAdvancedOptions) {
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
			shortName="Walnut";
			name="Walnut: PECAN-based Peptide Centric Analysis";
			image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/pecan_icon.png"));
		} else if (ProgramType.CASiL==program) {
			shortName="Thesaurus";
			name="Thesaurus: Phosphopeptide Positional Isomer Search Engine";
			image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/thesaurus_icon.png"));
		} else if (ProgramType.XCorDIA==program) {
			shortName="XCorDIA";
			name="XCorDIA: Peptide Searching with Cross Correlation";
			image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/mike_rotate_icon.png"));
		} else {
			if (isOffending==1) {
				shortName="MaiziepeDIA";
				name="Lindsay's MaiziepeDIA: Peptide Searching for DIA";
				image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/mazie_icon.png"));
			} else if (isOffending==2) {
				shortName="ChocopeDIA";
				name="MoMo's ChocopeDIA: Peptide Searching for DIA";
				image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/chocolate2.png"));
			} else if (isOffending==3) {
				shortName="StupiDIA";
				name="Austin's StupiDIA: Peptide Searching for DIA";
				image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/austin_icon.png"));
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
			try {
				// Replace: import com.apple.eawt.Application
				String className="com.apple.eawt.Application";
				Class<?> cls=Class.forName(className);

				// Replace: Application application=Application.getApplication();
				
				Object application=cls.newInstance().getClass().getMethod("getApplication").invoke(null);
				

				// Replace: application.setDockIconImage(image);
				application.getClass().getMethod("setDockIconImage", java.awt.Image.class).invoke(application, image.getImage());
				
			} catch (Exception e) {
				Logger.errorLine("Error setting Mac-specific properties");
				Logger.errorException(e);
			}
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
		f.setJMenuBar(panel.createMenus(program, enableAdvancedOptions));

		f.pack();
		
		Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		f.setSize(new Dimension(Math.max(600, Math.min(screenSize.width-20, 1250)), Math.max(600, Math.min(screenSize.height-20, 950))));
		f.setVisible(true);

		Runtime instance=Runtime.getRuntime();
		long mbOfMemory=(instance.maxMemory()/MemoryMonitor.mb);
		if (mbOfMemory<1000) {
			JOptionPane.showMessageDialog(f, "Warning, you only have "+mbOfMemory+" MB of memory allocated.\nPlease make sure you are running 64-bit Java!", "Warning, Low Memory!", JOptionPane.WARNING_MESSAGE, image);
		}
		

		if (!Main.isJavaVersionOK()) {
			int javaVersion=Main.getJavaVersion();
			String text=javaVersion<8?"lower":"higher";
			JOptionPane.showMessageDialog(f, "Warning, Java version ("+javaVersion+") is "+text+" than expected (8-16), execution may be unstable!", "Warning, Incorrect Java Version!", JOptionPane.WARNING_MESSAGE, image);
		}

		Logger.logLine(shortName+" Graphical Interface (version "+program.getVersion()+")");
		
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
		
		VersioningDetector.checkVersionGUI(program, f);
		
		return f;
	}
}