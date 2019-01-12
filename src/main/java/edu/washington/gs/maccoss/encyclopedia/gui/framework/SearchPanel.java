package edu.washington.gs.maccoss.encyclopedia.gui.framework;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.FileDialog;
import java.awt.FlowLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Optional;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;
import javax.swing.table.TableColumn;

import edu.washington.gs.maccoss.encyclopedia.Encyclopedia;
import edu.washington.gs.maccoss.encyclopedia.Pecanpie;
import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.CASiLSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCordiaSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.gui.dia.DIABrowserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.dia.FeatureGrapher;
import edu.washington.gs.maccoss.encyclopedia.gui.dia.LocalizationResultsBrowserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.dia.MultiResultsBrowserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.dia.PeptideExtractingBrowserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.dia.ResultsBrowserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.dia.WindowingSchemeWizard;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.library.AustinsSpecialEncyclopediaPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.library.CASiLParametersPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.library.EncyclopediaParametersPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.library.LindsaysSpecialEncyclopediaPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.library.MoMosSpecialEncyclopediaPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.pecan.PecanParametersPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.xcordia.XCorDIAParametersPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.AboutDialog;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessorTableModel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LogConsole;
import edu.washington.gs.maccoss.encyclopedia.gui.general.MemoryMonitor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.ProgressRenderer;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.Networking;

public class SearchPanel extends JPanel {
	private static final long serialVersionUID=1L;

	private static final ImageIcon openIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/fileopen.png"));
	private static final ImageIcon skylineIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/skyline_icon.png"));
	private static final ImageIcon openDBIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/filedb.png"));
	private static final ImageIcon convertDBIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/convertdb.png"));
	private static final ImageIcon libraryBrowserIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/encyclopedia_small_icon.png"));
	private static final ImageIcon diaBrowserIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/orbi_icon.png"));
	private static final ImageIcon peptideBrowserIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/peptide_icon.png"));
	private static final ImageIcon featureBrowserIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/feature_icon.png"));
	private static final ImageIcon helpIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/help_icon.png"));
	private static final ImageIcon windowSchemeIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/window_scheme_icon.png"));
	
	JobProcessorTableModel processorTableModel=new JobProcessorTableModel();
	
	private final JTabbedPane optionsTabs;
	//private final JCheckBox alignBetweenFiles;
	
	public SearchPanel(ProgramType program) {
		super(new BorderLayout());
	    
		setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
		setBackground(Color.white);

		JSplitPane split=new JSplitPane();

		optionsTabs=new JTabbedPane();
		
		if (ProgramType.Global==program||ProgramType.EncyclopeDIA==program) {
			try {
				EncyclopediaParametersPanel encyclopedia;
				switch (Networking.isOffendingAddress()) {
					case 1:
						encyclopedia=new LindsaysSpecialEncyclopediaPanel(this);
						break;
					case 2:
						encyclopedia=new MoMosSpecialEncyclopediaPanel(this);
						break;
					case 3:
						encyclopedia=new AustinsSpecialEncyclopediaPanel(this);
						break;
					default:
						encyclopedia=new EncyclopediaParametersPanel(this);
						break;
				}
				HashMap<String, String> map=SearchParameters.readPreferences();
				encyclopedia.setParameters(SearchParameterParser.parseParameters(map), map.get(Encyclopedia.TARGET_LIBRARY_TAG), map.get(Encyclopedia.BACKGROUND_FASTA_TAG));
				optionsTabs.addTab(encyclopedia.getProgramName(), encyclopedia.getSmallImage(), encyclopedia, encyclopedia.getProgramShortDescription());
			} catch (Exception e) {
				Logger.errorLine("Unexpected error reading saved parameters; using default parameters.");
				Logger.errorException(e);
			}
		}
		if (ProgramType.Global==program||ProgramType.CASiL==program) {
			try {
				CASiLParametersPanel CASiL=new CASiLParametersPanel(this);
				HashMap<String, String> map=CASiLSearchParameters.readPreferences();
				CASiLSearchParameters xcordiaParameters=CASiLSearchParameters.convertFromEncyclopeDIA(SearchParameterParser.parseParameters(map));
				CASiL.setParameters(xcordiaParameters, map.get(Encyclopedia.TARGET_LIBRARY_TAG), map.get(Encyclopedia.BACKGROUND_FASTA_TAG));
				optionsTabs.addTab(CASiL.getProgram().toString(), CASiL.getSmallImage(), CASiL, CASiL.getProgramShortDescription());
			} catch (Exception e) {
				Logger.errorLine("Unexpected error reading saved parameters; using default parameters.");
				Logger.errorException(e);
			}
		}
		if (ProgramType.Global==program||ProgramType.PecanPie==program||ProgramType.EncyclopeDIA==program) {
			PecanParametersPanel pecan=new PecanParametersPanel(this);
			try {
				HashMap<String, String> map=PecanSearchParameters.readPreferences();
				PecanSearchParameters parseParameters=PecanParameterParser.parseParameters(map);
				pecan.setParameters(parseParameters, map.get(Pecanpie.BACKGROUND_FASTA_TAG), map.get(Pecanpie.TARGET_FASTA_TAG));
				
			} catch (Exception e) {
				Logger.errorLine("Unexpected error reading saved parameters; using default parameters.");
				Logger.errorException(e);
			}
			optionsTabs.addTab(pecan.getProgram().toString(), pecan.getSmallImage(), pecan, pecan.getProgramShortDescription());
		}
		if (ProgramType.Global==program||ProgramType.XCorDIA==program) {
			XCorDIAParametersPanel xcordia=new XCorDIAParametersPanel(this);
			try {
				HashMap<String, String> map=XCordiaSearchParameters.readPreferences();
				XCordiaSearchParameters xcordiaParameters=XCordiaSearchParameters.convertFromPecan(PecanParameterParser.parseParameters(map));
				xcordia.setParameters(xcordiaParameters, map.get(Pecanpie.BACKGROUND_FASTA_TAG), map.get(Pecanpie.TARGET_FASTA_TAG));
				
			} catch (Exception e) {
				Logger.errorLine("Unexpected error reading saved parameters; using default parameters.");
				Logger.errorException(e);
			}
			optionsTabs.addTab(xcordia.getProgram().toString(), xcordia.getSmallImage(), xcordia, xcordia.getProgramShortDescription());
		}

		LogConsole console=new LogConsole();
		console.errorLine("Console:");
		Logger.addRecorder(console);

		MemoryMonitor memory=new MemoryMonitor();
		memory.start();
		
		JPanel optionsWrapper=new JPanel(new BorderLayout());
		optionsWrapper.setOpaque(true);
		optionsWrapper.setBackground(Color.white);
		optionsWrapper.add((Component)optionsTabs, BorderLayout.NORTH);
		optionsWrapper.add(console, BorderLayout.CENTER);
		optionsWrapper.add(memory, BorderLayout.SOUTH);
		
		split.setLeftComponent(optionsWrapper);

		JPanel files=new JPanel(new BorderLayout());
		JButton chooseFile=new JButton("Add MZML", openIcon);
		chooseFile.setToolTipText("Add MZML to analysis stack and process using current settings.");
		chooseFile.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				addMZML();
			}
		});
		
		//alignBetweenFiles=new JCheckBox("RT Align", true);
		//alignBetweenFiles.setToolTipText("Align retention times between files. Only uncheck for generating searchable chromatogram libraries where fractions don't share peptides.");
		
		JButton saveBlib=new JButton("Save BLIB", skylineIcon);
		saveBlib.setToolTipText("Save Skyline BLIB library.");
		saveBlib.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				saveBLIB(true);
			}
		});
		
		JButton saveChromElib=new JButton("Save Chromatogram Library", openDBIcon);
		saveChromElib.setToolTipText("Save chromatogram library as ELIB.");
		saveChromElib.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				saveELIB(false);
			}
		});
		
		JButton saveElib=new JButton("Save Quant Reports", openDBIcon);
		saveElib.setToolTipText("Save quantitative reports as ELIB.");
		saveElib.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				saveELIB(true);
			}
		});
		
		
		JPanel buttonPanel=new JPanel(new FlowLayout());
		buttonPanel.add(chooseFile);
		//buttonPanel.add(alignBetweenFiles);

		if (ProgramType.PecanPie!=program) {
			buttonPanel.add(saveChromElib);
			buttonPanel.add(saveElib);
		}
		
		buttonPanel.add(saveBlib);
		
		files.add(new JLabel("<html><p style=\"font-size:12px; font-family: Helvetica, sans-serif\"><b>Jobs: "), BorderLayout.WEST);
		files.add(buttonPanel, BorderLayout.EAST);
		
		JPanel filesWrapper=new JPanel(new BorderLayout());
		filesWrapper.setOpaque(true);
		filesWrapper.setBackground(Color.white);
		filesWrapper.add(files, BorderLayout.NORTH);
		
		JTable table=new JTable(processorTableModel);
        TableColumn column = table.getColumnModel().getColumn(1);
        column.setCellRenderer(new ProgressRenderer());
		filesWrapper.add(new JScrollPane(table), BorderLayout.CENTER);
		
		split.setRightComponent(filesWrapper);
		split.setDividerLocation(600);

		this.add(split, BorderLayout.CENTER);
	}
	
	public JMenuBar createMenus(ProgramType program) {
		JMenuBar bar=new JMenuBar();
		JMenu fileMenu=new JMenu("File");
		fileMenu.setMnemonic(KeyEvent.VK_F);
		bar.add(fileMenu);

		JMenuItem loadTarget=new JMenuItem("Load Target File", openIcon);
		loadTarget.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				loadTargetFile();
			}
		});
		loadTarget.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_L, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		fileMenu.add(loadTarget);

		JMenuItem openMZML=new JMenuItem("Open RAW File", openIcon);
		openMZML.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				addMZML();
			}
		});
		openMZML.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_R, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		fileMenu.add(openMZML);
		
		fileMenu.addSeparator();

		JMenuItem saveELIB=new JMenuItem("Save Quant Reports ELIB", openDBIcon);
		saveELIB.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				saveELIB(true);
			}
		});
		
		saveELIB.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));

		if (ProgramType.PecanPie!=program) {
			fileMenu.add(saveELIB);
		}

		JMenuItem saveBLIB=new JMenuItem("Save Skyline BLIB", skylineIcon);
		saveBLIB.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				saveBLIB(true);
			}
		});
		saveBLIB.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_B, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		fileMenu.add(saveBLIB);

		JMenu viewMenu=new JMenu("View");
		viewMenu.setMnemonic(KeyEvent.VK_V);

		if (ProgramType.PecanPie!=program) {
			bar.add(viewMenu);
		}

		JMenuItem launchBrowser=new JMenuItem("Launch ELIB Browser", libraryBrowserIcon);
		launchBrowser.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				launchElibBrowser();
			}
		});
		launchBrowser.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		viewMenu.add(launchBrowser);

		JMenuItem launchMultiBrowser=new JMenuItem("Launch Multi-ELIB Browser", libraryBrowserIcon);
		launchMultiBrowser.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				launchMultiElibBrowser();
			}
		});
		launchMultiBrowser.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_M, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));

		if (ProgramType.CASiL!=program) {
			viewMenu.add(launchMultiBrowser);
		}
		
		viewMenu.addSeparator();

		JMenuItem launchDIABrowser=new JMenuItem("Launch RAW File Browser", diaBrowserIcon);
		launchDIABrowser.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				launchDIABrowser();
			}
		});
		viewMenu.add(launchDIABrowser);

		JMenuItem launchPeptideBrowser=new JMenuItem("Launch Peptide Browser", peptideBrowserIcon);
		launchPeptideBrowser.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				launchPeptideBrowser();
			}
		});
		if (ProgramType.CASiL!=program) {
			viewMenu.add(launchPeptideBrowser);
		}

		JMenuItem launchFeatureBrowser=new JMenuItem("Launch Feature Browser", featureBrowserIcon);
		launchFeatureBrowser.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				launchFeatureBrowser();
			}
		});
		if (ProgramType.CASiL!=program) {
			viewMenu.add(launchFeatureBrowser);
		}

		JMenu convertMenu=new JMenu("Convert");
		convertMenu.setMnemonic(KeyEvent.VK_C);

		if (ProgramType.PecanPie!=program) {
			bar.add(convertMenu);
		}
		
		JMenuItem convertBLIB=new JMenuItem("Convert BLIB to Library", convertDBIcon);
		convertBLIB.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.convertBLIB(SearchPanel.this, getVisibleTab().getParameters());
			}
		});
		convertMenu.add(convertBLIB);

		JMenuItem convertMSP=new JMenuItem("Convert SPTXT/MSP to Library", convertDBIcon);
		convertMSP.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.convertMSP(SearchPanel.this, getVisibleTab().getParameters());
			}
		});
		convertMenu.add(convertMSP);

		JMenuItem convertTraML=new JMenuItem("Convert TraML to Library", convertDBIcon);
		convertTraML.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.convertTRAML(SearchPanel.this, getVisibleTab().getParameters());
			}
		});
		convertMenu.add(convertTraML);

		// FIXME NOT QUITE WORKING YET
//		JMenuItem convertOStsv=new JMenuItem("Convert OpenSWATH tsv to Library", convertDBIcon);
//		convertOStsv.addActionListener(new ActionListener() {
//			@Override
//			public void actionPerformed(ActionEvent e) {
//				SearchPanelUtilities.convertOpenSwathToELIB(SearchPanel.this, getVisibleTab().getParameters());
//			}
//		});
//		convertMenu.add(convertOStsv);
		
		convertMenu.addSeparator();
		
		JMenuItem convertELIBtoBLIB=new JMenuItem("Convert Library to BLIB", convertDBIcon);
		convertELIBtoBLIB.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.convertELIBtoBLIB(SearchPanel.this);
			}
		});
		convertMenu.add(convertELIBtoBLIB);
		
		JMenuItem convertELIBtoOpenSWATH=new JMenuItem("Convert Library to OpenSWATH tsv", convertDBIcon);
		convertELIBtoOpenSWATH.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.convertELIBtoOpenSWATH(SearchPanel.this, getVisibleTab().getParameters());
			}
		});
		convertMenu.add(convertELIBtoOpenSWATH);
		
		convertMenu.addSeparator();
		
		JMenuItem combineELIB=new JMenuItem("Combine Multiple Libraries", convertDBIcon);
		combineELIB.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.combineELIBs(SearchPanel.this);
			}
		});
		convertMenu.add(combineELIB);
		
		JMenuItem subsetELIB=new JMenuItem("Create Subset Library", convertDBIcon);
		subsetELIB.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.subsetELIB(SearchPanel.this);
			}
		});
		convertMenu.add(subsetELIB);
		
		JMenu helpMenu=new JMenu("Help");
		helpMenu.setMnemonic(KeyEvent.VK_H);
		bar.add(helpMenu);

		JMenuItem aboutMenuItem=new JMenuItem("About", helpIcon);
		aboutMenuItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				about();
			}
		});
		helpMenu.add(aboutMenuItem);

		helpMenu.addSeparator();

		JMenuItem windowSchemeItem=new JMenuItem("Window Scheme Wizard", windowSchemeIcon);
		windowSchemeItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				launchWindowingSchemeWizard();
			}
		});
		helpMenu.add(windowSchemeItem);
		
		return bar;
	}
	
	public void about() {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(SearchPanel.this);
		ParametersPanelInterface panel=getVisibleTab();
		AboutDialog.showAbout(frame, panel.getProgram(), panel.getImage());
	}
	
	public Collection<ParametersPanelInterface> getAllTabs() {
		ArrayList<ParametersPanelInterface> list=new ArrayList<ParametersPanelInterface>();
		for (int i=0; i<optionsTabs.getTabCount(); i++) {
			list.add((ParametersPanelInterface)optionsTabs.getComponentAt(i));
		}
		return list;
	}
	
	public ParametersPanelInterface getVisibleTab() {
		return (ParametersPanelInterface)optionsTabs.getSelectedComponent();
	}
	
	public void launchFeatureBrowser() {
		File featureFile=FileChooserPanel.getFiles(null, "Feature text files", new SimpleFilenameFilter("features.txt"), (JFrame)null, true)[0];

		if (featureFile!=null&&featureFile.exists()) {
			final JFrame dialog=new JFrame("Global Feature Browser");

			dialog.getContentPane().add(FeatureGrapher.graphFeatures(featureFile), BorderLayout.CENTER);
			dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
			dialog.pack(); 
			dialog.setSize(1900, 1030);
			dialog.setVisible(true);
		}
	}
	
	public void launchPeptideBrowser() {
		final JFrame dialog=new JFrame("Peptide/DIA Detection Browser");

		dialog.getContentPane().add(new PeptideExtractingBrowserPanel(getVisibleTab().getParameters()), BorderLayout.CENTER);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		dialog.pack(); 
		dialog.setSize(1900, 1030);
		dialog.setVisible(true);
	}
	
	public void launchElibBrowser() {
		if (getVisibleTab() instanceof CASiLParametersPanel) {
			launchLocalizationBrowser();
			return;
		}
		
		final JFrame dialog=new JFrame("ELIB/DIA Detection Browser");

		JMenuBar bar=new JMenuBar();
		JMenu fileMenu=new JMenu("File");
		fileMenu.setMnemonic(KeyEvent.VK_F);
		bar.add(fileMenu);
		
		final ResultsBrowserPanel browser=new ResultsBrowserPanel(getVisibleTab().getParameters());
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
		dialog.setJMenuBar(bar);
		
		dialog.getContentPane().add(browser, BorderLayout.CENTER);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		dialog.pack(); 
		dialog.setSize(1900, 1030);
		dialog.setVisible(true);
	}
	
	public void launchLocalizationBrowser() {
		final JFrame dialog=new JFrame("Localization Browser");

		JMenuBar bar=new JMenuBar();
		JMenu fileMenu=new JMenu("File");
		fileMenu.setMnemonic(KeyEvent.VK_F);
		bar.add(fileMenu);
		
		final LocalizationResultsBrowserPanel browser=new LocalizationResultsBrowserPanel(getVisibleTab().getParameters());
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
		dialog.setJMenuBar(bar);
		
		dialog.getContentPane().add(browser, BorderLayout.CENTER);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		dialog.pack(); 
		dialog.setSize(1900, 1030);
		dialog.setVisible(true);
	}
	
	public void launchWindowingSchemeWizard() {
		final JFrame dialog=new JFrame("Window Scheme Wizard");
		
		dialog.getContentPane().add(new WindowingSchemeWizard(), BorderLayout.CENTER);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		dialog.pack(); 
		dialog.setSize(1900, 1030);
		dialog.setVisible(true);
	}
	
	public void launchMultiElibBrowser() {
		final JFrame dialog=new JFrame("Multi ELIB/DIA Detection Browser");

		JMenuBar bar=new JMenuBar();
		JMenu fileMenu=new JMenu("File");
		fileMenu.setMnemonic(KeyEvent.VK_F);
		bar.add(fileMenu);
		
		final MultiResultsBrowserPanel browser=new MultiResultsBrowserPanel(getVisibleTab().getParameters());
		JMenuItem openElib=new JMenuItem("Open Multi ELIB...", libraryBrowserIcon);
		openElib.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				browser.askForLibrary();
			}
		});
		openElib.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		fileMenu.add(openElib);
		
		dialog.setJMenuBar(bar);
		
		dialog.getContentPane().add(browser, BorderLayout.CENTER);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		dialog.pack(); 
		dialog.setSize(1900, 1030);
		dialog.setVisible(true);
	}
	
	public void launchDIABrowser() {
		final JFrame dialog=new JFrame("RAW File Browser");

		JMenuBar bar=new JMenuBar();
		JMenu fileMenu=new JMenu("File");
		fileMenu.setMnemonic(KeyEvent.VK_F);
		bar.add(fileMenu);
		
		final DIABrowserPanel browser=new DIABrowserPanel(getVisibleTab().getParameters());
		
		JMenuItem rawBrowser=new JMenuItem("Open RAW File...", diaBrowserIcon);
		rawBrowser.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				browser.askForRaw();
			}
		});
		rawBrowser.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_R, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		fileMenu.add(rawBrowser);
		dialog.setJMenuBar(bar);
		
		dialog.getContentPane().add(browser, BorderLayout.CENTER);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		dialog.pack(); 
		dialog.setSize(1900, 1030);
		dialog.setVisible(true);
	}

	public void saveBLIB(boolean alignBetweenFiles) {
		JFrame frame = (JFrame)SwingUtilities.getRoot(SearchPanel.this);

		Optional<String> maybeError=getVisibleTab().canLoadData();
		if (maybeError.isPresent()) {
			JOptionPane.showMessageDialog(frame, maybeError.get());
		} else if (processorTableModel.getRowCount()==0) {
			JOptionPane.showMessageDialog(frame, "Please queue some RAW files first!");
			
		} else {
			FileDialog dialog=new FileDialog(frame, "Save a BLIB file", FileDialog.SAVE);
			dialog.setFilenameFilter(new SimpleFilenameFilter(".blib"));
			dialog.setVisible(true);
			if (dialog.getFiles()!=null&&dialog.getFiles().length>0) {
				File blibFile=dialog.getFiles()[0];
				String fileName=blibFile.getName();
				if (!fileName.toLowerCase().endsWith(".blib")) {
					blibFile=new File(blibFile.getParentFile(), fileName+".blib");

					if (blibFile.exists()) {
						// TODO ask if you want to overwrite this
						// updated file location!
					}
				}

				Logger.logLine("Adding BLIB export to queue for ["+blibFile.getAbsolutePath()+"]");
				SearchToBLIBJob job=new SearchToBLIBJob(blibFile, alignBetweenFiles, processorTableModel);
				if (job!=null) {
					processorTableModel.addJob(job);
				}
			}
		}
	}

//	public boolean isAlignedBetweenFiles() {
//		return alignBetweenFiles.isSelected();
//	}

	public void saveELIB(boolean alignBetweenFiles) {
		JFrame frame = (JFrame)SwingUtilities.getRoot(SearchPanel.this);

		Optional<String> maybeError=getVisibleTab().canLoadData();
		if (maybeError.isPresent()) {
			JOptionPane.showMessageDialog(frame, maybeError.get());
		} else if (processorTableModel.getRowCount()==0) {
			JOptionPane.showMessageDialog(frame, "Please queue some RAW files first!");
			
		} else {
			FileDialog dialog=new FileDialog(frame, "Save a ELIB file", FileDialog.SAVE);
			dialog.setFilenameFilter(new SimpleFilenameFilter(LibraryFile.ELIB));
			dialog.setVisible(true);
			if (dialog.getFiles()!=null&&dialog.getFiles().length>0) {
				File elibFile=dialog.getFiles()[0];
				String fileName=elibFile.getName();
				if (!fileName.toLowerCase().endsWith(LibraryFile.ELIB)) {
					elibFile=new File(elibFile.getParentFile(), fileName+LibraryFile.ELIB);

					if (elibFile.exists()) {
						// TODO ask if you want to overwrite this
						// updated file location!
					}
				}

				Logger.logLine("Adding ELIB export to queue for ["+elibFile.getAbsolutePath()+"]");
				SearchToELIBJob job=new SearchToELIBJob(elibFile, alignBetweenFiles, processorTableModel);
				if (job!=null) {
					processorTableModel.addJob(job);
				}
			}
		}
	}
	
	public void loadTargetFile() {
		getVisibleTab().askForSetupFile();
	}

	public void addMZML() {
		JFrame frame = (JFrame)SwingUtilities.getRoot(SearchPanel.this);
		
		Optional<String> maybeError=getVisibleTab().canLoadData();
		if (maybeError.isPresent()) {
			JOptionPane.showMessageDialog(frame, maybeError.get());
		} else {
			FileDialog dialog=new FileDialog(frame, "Select a RAW file", FileDialog.LOAD);
			dialog.setMultipleMode(true);
			dialog.setFilenameFilter(StripeFileGenerator.getFilenameFilter());
			dialog.setVisible(true);
			if (dialog.getFiles()!=null) {
				for (File file : dialog.getFiles()) {
					Logger.logLine("Adding mzML import to queue for ["+file.getAbsolutePath()+"]");
					getVisibleTab().getJob(file, processorTableModel);
				}
			}
		}
	}
}
