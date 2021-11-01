package edu.washington.gs.maccoss.encyclopedia;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.HashMap;
import java.util.Optional;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.dia.ResultsBrowserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.MemoryMonitor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingWorkerProgress;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;

/**
 * CASiL or CAPSLoc (Chromatogram Aligned Phospho Site Localizer)
 * @author searleb
 *
 */
public class PhosphoSiteLocalizerGUI {
	private static final String shortName="Thesaurus";
	private static final String name="Thesaurus: Chromatogram Aligned Site Localizing Search Engine";
	private static final ImageIcon image=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/mazie_icon.png"));
	private static final ImageIcon libraryBrowserIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/encyclopedia_small_icon.png"));
	private static final ImageIcon diaBrowserIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/orbi_icon.png"));
	private static final ImageIcon skylineIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/skyline_icon.png"));

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
		
		final ResultsBrowserPanel browser=new ResultsBrowserPanel(getParameters());
		

		JMenuBar bar=new JMenuBar();
		JMenu fileMenu=new JMenu("File");
		fileMenu.setMnemonic(KeyEvent.VK_F);
		bar.add(fileMenu);

		JMenuItem convertBLIB=new JMenuItem("Convert BLIB to ELIB", skylineIcon);
		convertBLIB.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				convertBLIB(f);
			}
		});
		convertBLIB.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_B, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		fileMenu.add(convertBLIB);
		
		JMenuItem openElib=new JMenuItem("Open ELIB...", libraryBrowserIcon);
		openElib.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				browser.askForLibrary();
			}
		});
		openElib.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		fileMenu.add(openElib);
		
		JMenuItem rawBrowser=new JMenuItem("Open RAW File...", diaBrowserIcon);
		rawBrowser.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				browser.askForRaw();
			}
		});
		rawBrowser.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_R, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		fileMenu.add(rawBrowser);
		f.setJMenuBar(bar);

		f.getContentPane().add(browser, BorderLayout.CENTER);

		f.pack();
		f.setSize(new Dimension(1250, 750));
		f.setVisible(true);

		Runtime instance=Runtime.getRuntime();
		long mbOfMemory=(instance.maxMemory()/MemoryMonitor.mb);
		if (mbOfMemory<1000) {
			JOptionPane.showMessageDialog(f, "Warning, you only have "+mbOfMemory+" MB of memory allocated.\nPlease make sure you are running 64-bit Java!", "Warning, Low Memory!", JOptionPane.WARNING_MESSAGE, image);
		}

		Logger.logLine(shortName+" Graphical Interface");
	}
	
	public static void convertBLIB(JFrame frame) {
		final JDialog dialog=new JDialog(frame, "Convert BLIB to ELIB", true);
		
		final FileChooserPanel blibFileChooser=new FileChooserPanel(null, "BLIB", new SimpleFilenameFilter(".blib"), true);
		final FileChooserPanel fastaFileChooser=new FileChooserPanel(null, "FASTA", new SimpleFilenameFilter(".fas", ".fasta"), true);

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(blibFileChooser);
		options.add(fastaFileChooser);
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				dialog.setVisible(false);
				dialog.dispose();

				final File blibFile=blibFileChooser.getFile();
				final File fastaFile=fastaFileChooser.getFile();
				
				if (blibFile!=null&&blibFile.exists()&&fastaFile!=null&&fastaFile.exists()) {
					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame)SwingUtilities.getWindowAncestor(frame), "Please wait...", "Reading BLIB File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							BlibToLibraryConverter.convert(blibFile, Optional.ofNullable(null), fastaFile, true, getParameters());
							Logger.logLine("Finished reading "+blibFile.getName());
							return Nothing.NOTHING;
						}
						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
					
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify a BLIB and a FASTA file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, skylineIcon);
				}
			}
		});
		buttons.add(okButton);
		JButton cancelButton=new JButton("Cancel");
		cancelButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				dialog.setVisible(false);
				dialog.dispose();
			}
		});
		buttons.add(cancelButton);
		
		JPanel mainpane=new JPanel(new BorderLayout());
		mainpane.add(options, BorderLayout.CENTER);
		mainpane.add(buttons, BorderLayout.SOUTH);
		mainpane.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10), BorderFactory.createTitledBorder("Parameters:")));
		
		dialog.getContentPane().add(mainpane, BorderLayout.CENTER);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		dialog.pack(); 
		dialog.setSize(500, 200);
		dialog.setVisible(true);
	}
	
	public static SearchParameters getParameters() {
		HashMap<String, String> map=SearchParameterParser.getDefaultParameters();
		map.put("-deconvoluteOverlappingWindows", "true");
		map.put("-fixed", "");
		map.put("-ptol", "50");
		map.put("-ftol", "50");
		map.put("-lftol", "50");
		SearchParameters parameters=SearchParameterParser.parseParameters(map);
		return parameters;
	}
}
