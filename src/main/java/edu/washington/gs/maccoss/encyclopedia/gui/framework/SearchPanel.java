package edu.washington.gs.maccoss.encyclopedia.gui.framework;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.FileDialog;
import java.awt.FlowLayout;
import java.awt.Insets;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Optional;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.KeyStroke;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.table.TableColumn;

import edu.washington.gs.maccoss.encyclopedia.Encyclopedia;
import edu.washington.gs.maccoss.encyclopedia.EncyclopediaTwo;
import edu.washington.gs.maccoss.encyclopedia.Pecanpie;
import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ThesaurusSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.scribe.ScribeSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCordiaSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.parameters.InstrumentSpecificSearchParameters;
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
import edu.washington.gs.maccoss.encyclopedia.gui.dia.interactive.ChromatogrindrPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.library.AustinsSpecialEncyclopediaPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.library.EncyclopediaParametersPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.library.EncyclopediaTwoParametersPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.library.LindsaysSpecialEncyclopediaPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.library.MoMosSpecialEncyclopediaPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.library.ThesaurusParametersPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.pecan.PecanParametersPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.scribe.ScribeParametersPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.xcordia.XCorDIAParametersPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.AboutDialog;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessorTableModel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LabeledComponent;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LogConsole;
import edu.washington.gs.maccoss.encyclopedia.gui.general.MemoryMonitor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.ProgressRenderer;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.jobs.JobProcessor;
import edu.washington.gs.maccoss.encyclopedia.jobs.SearchToBLIBJob;
import edu.washington.gs.maccoss.encyclopedia.jobs.SearchToELIBJob;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.Networking;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

public class SearchPanel extends JPanel {
	private static final long serialVersionUID=1L;

	private static final ImageIcon openIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/fileopen.png"));
	private static final ImageIcon skylineIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/skyline_icon.png"));
	private static final ImageIcon openDBIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/filedb.png"));
	private static final ImageIcon convertDBIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/convertdb.png"));
	private static final ImageIcon processingIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/processing_icon.png"));
	private static final ImageIcon libraryBrowserIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/library_small_icon.png"));
	private static final ImageIcon diaBrowserIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/orbi_icon.png"));
	private static final ImageIcon peptideBrowserIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/peptide_icon.png"));
	private static final ImageIcon featureBrowserIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/feature_icon.png"));
	private static final ImageIcon grindrIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/grindr_icon.png"));
	private static final ImageIcon helpIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/help_icon.png"));
	private static final ImageIcon windowSchemeIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/window_scheme_icon.png"));
	private static final ImageIcon xmlFileIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/xml_file_icon.png"));

	private static final int numberOfCores=Runtime.getRuntime().availableProcessors();
	private final SpinnerModel numberOfJobs=new SpinnerNumberModel(numberOfCores, 1, numberOfCores, 1);
	private final JComboBox<InstrumentSpecificSearchParameters> instrumentCombo=new JComboBox<InstrumentSpecificSearchParameters>(InstrumentSpecificSearchParameters.INSTRUMENTS);
	private final JComboBox<DigestionEnzyme> enzymeCombo=new JComboBox<DigestionEnzyme>(new Vector<>(DigestionEnzyme.getAvailableEnzymes()));
	
	JobProcessorTableModel processorTableModel=new JobProcessorTableModel();
	
	private final JTabbedPane engineSpecificParameters;
	//private final JCheckBox alignBetweenFiles;
	

	public void setParameters(SearchParameters params) {
		instrumentCombo.setSelectedItem(params.getInstrument());
		enzymeCombo.setSelectedItem(params.getEnzyme());
		numberOfJobs.setValue(params.getNumberOfThreadsUsed());
	}
	
	public SearchPanel(ProgramType program, boolean enableAdvancedOptions) {
		super(new BorderLayout());
	    
		setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
		setBackground(Color.white);

		JSplitPane split=new JSplitPane();

		JPanel generalParameters=new JPanel();
		generalParameters.setBackground(Color.white);
		generalParameters.setLayout(new BoxLayout(generalParameters, BoxLayout.PAGE_AXIS));
		
		generalParameters.add(new JSeparator(JSeparator.HORIZONTAL));
		generalParameters.add(new LabeledComponent("<p style=\"font-size:12px; font-family: Helvetica, sans-serif\"><b>General Parameters", new JLabel()));
		generalParameters.add(new LabeledComponent("Instrument", instrumentCombo));
		generalParameters.add(new LabeledComponent("Enzyme", enzymeCombo));
		generalParameters.add(new LabeledComponent("Number of Cores", new JSpinner(numberOfJobs)));
		
		// to make it have the same margins as the JTabbedPane
		Insets insets = UIManager.getInsets("TabbedPane.tabInsets");
		generalParameters.setBorder(BorderFactory.createEmptyBorder(insets.top, insets.left, insets.bottom, insets.right));

		// set up JTabbedPane
		engineSpecificParameters=new JTabbedPane();

		HashMap<String, String> map;
		SearchParameters params;
		try {
			map=SearchParameters.readPreferences();
			params=SearchParameterParser.parseParameters(map);
			setParameters(params);
		} catch (Exception e) {
			Logger.errorLine("Unexpected error reading saved parameters; using default parameters.");
			Logger.errorException(e);
			map=SearchParameterParser.getDefaultParameters();
			params=SearchParameterParser.parseParameters(map);
		}

		if ((ProgramType.Global==program||ProgramType.EncyclopeDIA==program)&&false) {
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
				encyclopedia.setParameters(params, map.get(Encyclopedia.TARGET_LIBRARY_TAG), map.get(Encyclopedia.BACKGROUND_FASTA_TAG));
				engineSpecificParameters.addTab(encyclopedia.getProgramName(), encyclopedia.getSmallImage(), encyclopedia, encyclopedia.getProgramShortDescription());
			} catch (Exception e) {
				Logger.errorLine("Unexpected error reading saved parameters; using default parameters.");
				Logger.errorException(e);
			}
		}
		if ((ProgramType.Global==program||ProgramType.EncyclopeDIA==program)) {
			try {
				EncyclopediaTwoParametersPanel encyclopedia=new EncyclopediaTwoParametersPanel(this);
				
				encyclopedia.setParameters(params, map.get(EncyclopediaTwo.PREALIGNMENT_LIBRARY_TAG), map.get(EncyclopediaTwo.TARGET_LIBRARY_TAG), map.get(EncyclopediaTwo.BACKGROUND_FASTA_TAG));
				engineSpecificParameters.addTab(encyclopedia.getProgramName(), encyclopedia.getSmallImage(), encyclopedia, encyclopedia.getProgramShortDescription());
			} catch (Exception e) {
				Logger.errorLine("Unexpected error reading saved parameters; using default parameters.");
				Logger.errorException(e);
			}
		}
		if (ProgramType.Global==program||ProgramType.CASiL==program||ProgramType.EncyclopeDIA==program) {
			try {
				ThesaurusParametersPanel CASiL=new ThesaurusParametersPanel(this);
				map=ThesaurusSearchParameters.readPreferences();
				ThesaurusSearchParameters thesaurusParameters=ThesaurusSearchParameters.parseParameters(map);
				CASiL.setParameters(thesaurusParameters, map.get(Encyclopedia.TARGET_LIBRARY_TAG), map.get(Encyclopedia.BACKGROUND_FASTA_TAG));
				engineSpecificParameters.addTab(CASiL.getProgram().toString(), CASiL.getSmallImage(), CASiL, CASiL.getProgramShortDescription());
			} catch (Exception e) {
				Logger.errorLine("Unexpected error reading saved parameters; using default parameters.");
				Logger.errorException(e);
			}
		}
		if (ProgramType.Global==program||ProgramType.PecanPie==program||ProgramType.EncyclopeDIA==program) {
			PecanParametersPanel pecan=new PecanParametersPanel(this);
			try {
				map=PecanSearchParameters.readPreferences();
				PecanSearchParameters parseParameters=PecanParameterParser.parseParameters(map);
				pecan.setParameters(parseParameters, map.get(Pecanpie.BACKGROUND_FASTA_TAG), map.get(Pecanpie.TARGET_FASTA_TAG));
				
			} catch (Exception e) {
				Logger.errorLine("Unexpected error reading saved parameters; using default parameters.");
				Logger.errorException(e);
			}
			engineSpecificParameters.addTab(pecan.getProgram().toString(), pecan.getSmallImage(), pecan, pecan.getProgramShortDescription());
		}
		if (ProgramType.Global==program||ProgramType.XCorDIA==program) {
			XCorDIAParametersPanel xcordia=new XCorDIAParametersPanel(this);
			try {
				map=XCordiaSearchParameters.readPreferences();
				XCordiaSearchParameters xcordiaParameters=XCordiaSearchParameters.convertFromPecan(PecanParameterParser.parseParameters(map));
				xcordia.setParameters(xcordiaParameters, map.get(Pecanpie.BACKGROUND_FASTA_TAG), map.get(Pecanpie.TARGET_FASTA_TAG));
				
			} catch (Exception e) {
				Logger.errorLine("Unexpected error reading saved parameters; using default parameters.");
				Logger.errorException(e);
			}
			engineSpecificParameters.addTab(xcordia.getProgram().toString(), xcordia.getSmallImage(), xcordia, xcordia.getProgramShortDescription());
		}
		if ((ProgramType.Global==program)||ProgramType.Scribe==program||ProgramType.EncyclopeDIA==program) {
			ScribeParametersPanel scribe=new ScribeParametersPanel(this);
			try {
				map=ScribeSearchParameters.readPreferences();
				ScribeSearchParameters scribeParameters=ScribeSearchParameters.convertFromEncyclopeDIA(SearchParameterParser.parseParameters(map));
				scribe.setParameters(scribeParameters, map.get(Encyclopedia.TARGET_LIBRARY_TAG), map.get(Encyclopedia.BACKGROUND_FASTA_TAG));
				
			} catch (Exception e) {
				Logger.errorLine("Unexpected error reading saved parameters; using default parameters.");
				Logger.errorException(e);
			}
			engineSpecificParameters.addTab(scribe.getProgram().toString(), scribe.getSmallImage(), scribe, scribe.getProgramShortDescription());
		}

		LogConsole console=new LogConsole();
		console.errorLine("Console:");
		Logger.addRecorder(console);

		MemoryMonitor memory=new MemoryMonitor();
		memory.start();
		
		JPanel settingsWrapper=new JPanel(new BorderLayout());
		settingsWrapper.setOpaque(true);
		settingsWrapper.setBackground(Color.white);
		settingsWrapper.add((Component)engineSpecificParameters, BorderLayout.CENTER);
		settingsWrapper.add(generalParameters, BorderLayout.SOUTH);
		
		JPanel optionsWrapper=new JPanel(new BorderLayout());
		optionsWrapper.setOpaque(true);
		optionsWrapper.setBackground(Color.white);
		optionsWrapper.add(settingsWrapper, BorderLayout.NORTH);
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

		if (ProgramType.PecanPie==program) {
			buttonPanel.add(saveBlib);
		} else if (ProgramType.Scribe==program) {
			buttonPanel.add(saveElib);
		} else {
			buttonPanel.add(saveChromElib);
			buttonPanel.add(saveElib);
			buttonPanel.add(saveBlib);
		}
		
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
	
	public JMenuBar createMenus(ProgramType program, boolean enableAdvancedOptions) {
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

		JMenuItem clearJobs=new JMenuItem("Clear Job List");
		clearJobs.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				clearJobs();
			}
		});
		
		clearJobs.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_X, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		fileMenu.add(clearJobs);
		
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

		JMenuItem chromatogrindrItem=new JMenuItem("Launch Chromatogrind Browser", grindrIcon);
		chromatogrindrItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				Logger.logLine("Launching Interactive Chromatogram Visualizer");
				ChromatogrindrPanel.launchBrowserPanel(new ChromatogrindrPanel());
			}
		});
		chromatogrindrItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_G, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		viewMenu.add(chromatogrindrItem);

		JMenu convertMenu=new JMenu("Convert");
		convertMenu.setMnemonic(KeyEvent.VK_C);

		if (ProgramType.PecanPie!=program) {
			bar.add(convertMenu);
		}
		
		JMenuItem fastaToProsit=new JMenuItem("Create Prosit Library from FASTA with Koina", convertDBIcon);
		fastaToProsit.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.convertFastaForProsit(SearchPanel.this, getJobProcessor(), getVisibleTab().getParameters());
			}
		});
		convertMenu.add(fastaToProsit);
		
		JMenuItem libraryToProsit=new JMenuItem("Create Prosit Library from Library with Koina", convertDBIcon);
		libraryToProsit.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.convertLibraryForProsit(SearchPanel.this, getJobProcessor(), getVisibleTab().getParameters());
			}
		});
		convertMenu.add(libraryToProsit);
		
		JMenuItem fastaToMS2PIP=new JMenuItem("Create MS2PIP PEPREC from FASTA", convertDBIcon);
		fastaToMS2PIP.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.convertFastaForMS2PIP(SearchPanel.this);
			}
		});
		convertMenu.add(fastaToMS2PIP);
		
		JMenuItem libraryToMS2PIP=new JMenuItem("Create MS2PIP PEPREC from Library", convertDBIcon);
		libraryToMS2PIP.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.convertLibraryForMS2PIP(SearchPanel.this);
			}
		});
		convertMenu.add(libraryToMS2PIP);
		
		convertMenu.addSeparator();
		
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

		JMenuItem convertSpectronaut=new JMenuItem("Convert Prosit/Spectronaut CSV/XLS to Library", convertDBIcon);
		convertSpectronaut.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.convertSpectronaut(SearchPanel.this, getVisibleTab().getParameters());
			}
		});
		convertMenu.add(convertSpectronaut);

		JMenuItem convertMaxquant=new JMenuItem("Convert Maxquant msms.txt to Library", convertDBIcon);
		convertMaxquant.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.convertMaxQuantMSMSTXT(SearchPanel.this, getVisibleTab().getParameters());
			}
		});
		convertMenu.add(convertMaxquant);

		JMenuItem convertTraML=new JMenuItem("Convert TraML to Library", convertDBIcon);
		convertTraML.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.convertTRAML(SearchPanel.this, getVisibleTab().getParameters());
			}
		});
		convertMenu.add(convertTraML);

		JMenuItem convertOStsv=new JMenuItem("Convert OpenSWATH tsv to Library", convertDBIcon);
		convertOStsv.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.convertOpenSwathToELIB(SearchPanel.this, getVisibleTab().getParameters());
			}
		});
		convertMenu.add(convertOStsv);

		JMenuItem convertMS2PIP=new JMenuItem("Convert MS2PIP csv to Library", convertDBIcon);
		convertMS2PIP.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.convertMS2PIPToELIB(SearchPanel.this, getVisibleTab().getParameters());
			}
		});
		convertMenu.add(convertMS2PIP);
		
		convertMenu.addSeparator();
		
		JMenuItem convertELIBtoBLIB=new JMenuItem("Convert Library to BLIB", convertDBIcon);
		convertELIBtoBLIB.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.convertELIBtoBLIB(SearchPanel.this);
			}
		});
		convertMenu.add(convertELIBtoBLIB);
		
		JMenuItem convertELIBtoMSP=new JMenuItem("Convert Library to NIST MSP", convertDBIcon);
		convertELIBtoMSP.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.convertELIBtoMSP(SearchPanel.this, getVisibleTab().getParameters());
			}
		});
		convertMenu.add(convertELIBtoMSP);
		
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
				SearchPanelUtilities.combineELIBs(SearchPanel.this, getVisibleTab().getParameters(), getJobProcessor());
			}
		});
		convertMenu.add(combineELIB);
		
		JMenuItem correctELIBRTs=new JMenuItem("Correct Library Retention Times", convertDBIcon);
		correctELIBRTs.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.correctLibraryRTs(SearchPanel.this, getVisibleTab().getParameters(), getJobProcessor());
			}
		});
		if (enableAdvancedOptions) {
			correctELIBRTs.setText("HIDDEN: "+correctELIBRTs.getText());
			convertMenu.add(correctELIBRTs);
		}
		
		JMenuItem subsetELIB=new JMenuItem("Create Subset Library", convertDBIcon);
		subsetELIB.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.subsetELIB(SearchPanel.this);
			}
		});
		convertMenu.add(subsetELIB);

		JMenu dataMenu=new JMenu("Data");
		dataMenu.setMnemonic(KeyEvent.VK_D);
		bar.add(dataMenu);

		JMenuItem saveDriverItem=new JMenuItem("Save XML driver file", xmlFileIcon);
		saveDriverItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.saveDriverFile(dataMenu, getVisibleTab().getParameters(), getJobProcessor());
			}
		});
		dataMenu.add(saveDriverItem);
		
		JMenuItem loadDriverItem=new JMenuItem("Load XML driver file", xmlFileIcon);
		loadDriverItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.loadDriverFile(dataMenu, getVisibleTab().getParameters(), getJobProcessor());
			}
		});
		dataMenu.add(loadDriverItem);
		
		JMenuItem mzmlPreprocessorItem=new JMenuItem("Preprocess mzMLs", convertDBIcon);
		mzmlPreprocessorItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.preprocessMZMLs(dataMenu, getVisibleTab().getParameters());
			}
		});
		dataMenu.add(mzmlPreprocessorItem);

		JMenuItem mzmlMergerItem=new JMenuItem("Combine Gas Phase Fractions", convertDBIcon);
		mzmlMergerItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.combineMZMLs(dataMenu, getJobProcessor(), getVisibleTab().getParameters());
			}
		});
		dataMenu.add(mzmlMergerItem);

		JMenuItem elibSeperatorItem=new JMenuItem("Extract Sample-Specific Libraries from ELIB", convertDBIcon);
		elibSeperatorItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.extractSampleSpecificDLIBs(dataMenu, getVisibleTab().getParameters(), getJobProcessor());
			}
		});
		if (enableAdvancedOptions) {
			elibSeperatorItem.setText("HIDDEN: "+elibSeperatorItem.getText());
			dataMenu.add(elibSeperatorItem);
		}
		
		JMenuItem subsetDIA=new JMenuItem("Create Subset mzML", convertDBIcon);
		subsetDIA.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.subsetDIA(SearchPanel.this, getVisibleTab().getParameters());
			}
		});
		if (enableAdvancedOptions) {
			subsetDIA.setText("HIDDEN: "+subsetDIA.getText());
			dataMenu.add(subsetDIA);
		}

		JMenuItem dilutonCurveFitterItem=new JMenuItem("Create PRM assay from DIA...", processingIcon);
		dilutonCurveFitterItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SearchPanelUtilities.launchDIAtoPRMDialog(SearchPanel.this);
			}
		});
		dataMenu.add(dilutonCurveFitterItem);

		JMenuItem toggleScoringSystemItem=new JMenuItem("Toggle EncyclopeDIA Scoring System", libraryBrowserIcon);
		toggleScoringSystemItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				EncyclopediaScoringFactory.USE_LEGACY_SCORING_SYSTEM=!EncyclopediaScoringFactory.USE_LEGACY_SCORING_SYSTEM;
				JOptionPane.showMessageDialog(SearchPanel.this,"USE_LEGACY_SCORING_SYSTEM set to "+EncyclopediaScoringFactory.USE_LEGACY_SCORING_SYSTEM,"EncyclopeDIA Scoring System", JOptionPane.PLAIN_MESSAGE);
			}
		});
		if (enableAdvancedOptions) {
			toggleScoringSystemItem.setText("HIDDEN: "+toggleScoringSystemItem.getText());
			dataMenu.add(toggleScoringSystemItem);
		}
		
		dataMenu.addSeparator();

		JMenuItem windowSchemeItem=new JMenuItem("Window Scheme Wizard", windowSchemeIcon);
		windowSchemeItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				launchWindowingSchemeWizard();
			}
		});
		dataMenu.add(windowSchemeItem);
		
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
		
		return bar;
	}
	
	public void about() {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(SearchPanel.this);
		ParametersPanelInterface panel=getVisibleTab();
		AboutDialog.showAbout(frame, panel.getProgram(), panel.getImage());
	}
	
	public Collection<ParametersPanelInterface> getAllTabs() {
		ArrayList<ParametersPanelInterface> list=new ArrayList<ParametersPanelInterface>();
		for (int i=0; i<engineSpecificParameters.getTabCount(); i++) {
			list.add((ParametersPanelInterface)engineSpecificParameters.getComponentAt(i));
		}
		return list;
	}
	
	public ParametersPanelInterface getVisibleTab() {
		return (ParametersPanelInterface)engineSpecificParameters.getSelectedComponent();
	}
	
	public void launchFeatureBrowser() {
		File[] featureFiles=FileChooserPanel.getFiles(null, "Feature text files", new SimpleFilenameFilter("features.txt", "features.txt.unsorted", ".pin"), (JFrame)null, true);

		if (featureFiles!=null&&featureFiles.length>0&&featureFiles[0].exists()) {
			final JFrame dialog=new JFrame("Global Feature Browser");

			dialog.getContentPane().add(FeatureGrapher.graphFeatures(featureFiles[0]), BorderLayout.CENTER);
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
		if (getVisibleTab() instanceof ThesaurusParametersPanel) {
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
				browser.askForResults();
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
	
	public JobProcessor getJobProcessor() {
		return processorTableModel;
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
					try {
						getVisibleTab().getJob(file, processorTableModel);
					} catch (Exception e) {
						JOptionPane.showMessageDialog(frame, e);
					}
				}
			}
		}
	}
	
	public void clearJobs() {
		processorTableModel.clearJobs();
	}
	
	public InstrumentSpecificSearchParameters getInstrument() {
		return (InstrumentSpecificSearchParameters)instrumentCombo.getSelectedItem();
	}
	
	public DigestionEnzyme getEnzyme() {
		return (DigestionEnzyme)enzymeCombo.getSelectedItem();
	}
	
	public int getNumberOfJobs() {
		return ((Integer)numberOfJobs.getValue()).intValue();
	}
}
