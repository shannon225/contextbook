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
import java.io.IOException;
import java.sql.SQLException;
import java.util.Optional;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
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

import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MSPReader;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.library.EncyclopediaParametersPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.pecan.PecanParametersPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessorTableModel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LogConsole;
import edu.washington.gs.maccoss.encyclopedia.gui.general.MemoryMonitor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.ProgressRenderer;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class SearchPanel extends JPanel {
	private static final long serialVersionUID=1L;

	private static final ImageIcon openIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/fileopen.png"));
	private static final ImageIcon skylineIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/skyline_icon.png"));
	private static final ImageIcon openDBIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/filedb.png"));
	private static final ImageIcon convertDBIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/convertdb.png"));
	
	JobProcessorTableModel processorTableModel=new JobProcessorTableModel();
	
	private final JTabbedPane optionsTabs;
	
	public SearchPanel(boolean pecanpie) {
		super(new BorderLayout());
	    
		setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
		setBackground(Color.white);

		JSplitPane split=new JSplitPane();

		optionsTabs=new JTabbedPane();
		if (!pecanpie) {
			optionsTabs.addTab("Encyclopedia", EncyclopediaParametersPanel.smallimage, new EncyclopediaParametersPanel(), "Encyclopedia Library Search");
		}
		optionsTabs.addTab("Pecan", PecanParametersPanel.smallimage, new PecanParametersPanel(), "Pecan Peptide Search");
		
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
		chooseFile.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				addMZML();
			}
		});
		
		JButton saveBlib=new JButton("Save BLIB", skylineIcon);
		saveBlib.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				saveBLIB();
			}
		});
		
		JButton saveElib=new JButton("Save Library", openDBIcon);
		saveElib.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				saveELIB();
			}
		});
		
		JPanel buttonPanel=new JPanel(new FlowLayout());
		buttonPanel.add(chooseFile);
		buttonPanel.add(saveElib);
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
	
	public JMenuBar createMenus() {
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

		JMenuItem openMZML=new JMenuItem("Open MZML", openIcon);
		openMZML.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				addMZML();
			}
		});
		openMZML.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		fileMenu.add(openMZML);

		JMenuItem saveELIB=new JMenuItem("Save Library", openDBIcon);
		saveELIB.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				saveELIB();
			}
		});
		saveELIB.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		fileMenu.add(saveELIB);

		JMenuItem saveBLIB=new JMenuItem("Save BLIB", skylineIcon);
		saveBLIB.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				saveBLIB();
			}
		});
		saveBLIB.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_B, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		fileMenu.add(saveBLIB);

		JMenu convertMenu=new JMenu("Convert");
		convertMenu.setMnemonic(KeyEvent.VK_C);
		bar.add(convertMenu);

		JMenuItem convertBLIB=new JMenuItem("Convert BLIB to Library", convertDBIcon);
		convertBLIB.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				convertBLIB();
			}
		});
		convertMenu.add(convertBLIB);

		JMenuItem convertMSP=new JMenuItem("Convert MSP to Library", convertDBIcon);
		convertMSP.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				convertMSP();
			}
		});
		convertMenu.add(convertMSP);
		
		return bar;
	}
	
	public ParametersPanelInterface getVisibleTab() {
		return (ParametersPanelInterface)optionsTabs.getSelectedComponent();
	}
	
	public void convertBLIB() {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(SearchPanel.this);
		final JDialog dialog=new JDialog(frame, "Convert BLIB to Library", true);
		
		final FileChooserPanel blibFileChooser=new FileChooserPanel(null, "BLIB", new SimpleFilenameFilter(".blib"), true);
		final FileChooserPanel iRTFileChooser=new FileChooserPanel(null, "IRT Database", new SimpleFilenameFilter(".irtdb"), false);
		final FileChooserPanel fastaFileChooser=new FileChooserPanel(null, "FASTA", new SimpleFilenameFilter(".fas", ".fasta"), true);

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(blibFileChooser);
		options.add(iRTFileChooser);
		options.add(fastaFileChooser);
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				dialog.setVisible(false);
				dialog.dispose();

				File blibFile=blibFileChooser.getFile();
				File irtFile=iRTFileChooser.getFile();
				File fastaFile=fastaFileChooser.getFile();
				
				if (blibFile!=null&&blibFile.exists()&&fastaFile!=null&&fastaFile.exists()) {
					BlibToLibraryConverter.convert(blibFile, Optional.ofNullable(irtFile), fastaFile);
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify a BLIB and a FASTA file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
	
	public void convertMSP() {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(SearchPanel.this);
		final JDialog dialog=new JDialog(frame, "Convert NIST MSP to Library", true);
		
		final FileChooserPanel blibFileChooser=new FileChooserPanel(null, "MSP", new SimpleFilenameFilter(".msp"), true);
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

				File blibFile=blibFileChooser.getFile();
				File fastaFile=fastaFileChooser.getFile();
				
				if (blibFile!=null&&blibFile.exists()&&fastaFile!=null&&fastaFile.exists()) {
					try {
						MSPReader.convertMSP(blibFile, fastaFile);
					} catch (IOException ioe) {
						throw new EncyclopediaException("ELIB writing IO error!", ioe);
					} catch (SQLException sqle) {
						throw new EncyclopediaException("ELIB writing SQL error!", sqle);
					}
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify a MSP and a FASTA file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
		dialog.setSize(500, 170);
		dialog.setVisible(true);
	}

	public void saveBLIB() {
		JFrame frame = (JFrame)SwingUtilities.getRoot(SearchPanel.this);

		Optional<String> maybeError=getVisibleTab().canLoadData();
		if (maybeError.isPresent()) {
			JOptionPane.showMessageDialog(frame, maybeError.get());
		} else if (processorTableModel.getRowCount()==0) {
			JOptionPane.showMessageDialog(frame, "Please queue some MZMLs first!");
			
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

				SearchToBLIBJob job=new SearchToBLIBJob(blibFile, processorTableModel);
				if (job!=null) {
					processorTableModel.addJob(job);
				}
			}
		}
	}

	public void saveELIB() {
		JFrame frame = (JFrame)SwingUtilities.getRoot(SearchPanel.this);

		Optional<String> maybeError=getVisibleTab().canLoadData();
		if (maybeError.isPresent()) {
			JOptionPane.showMessageDialog(frame, maybeError.get());
		} else if (processorTableModel.getRowCount()==0) {
			JOptionPane.showMessageDialog(frame, "Please queue some MZMLs first!");
			
		} else {
			FileDialog dialog=new FileDialog(frame, "Save a ELIB file", FileDialog.SAVE);
			dialog.setFilenameFilter(new SimpleFilenameFilter(".elib"));
			dialog.setVisible(true);
			if (dialog.getFiles()!=null&&dialog.getFiles().length>0) {
				File elibFile=dialog.getFiles()[0];
				String fileName=elibFile.getName();
				if (!fileName.toLowerCase().endsWith(".elib")) {
					elibFile=new File(elibFile.getParentFile(), fileName+".elib");

					if (elibFile.exists()) {
						// TODO ask if you want to overwrite this
						// updated file location!
					}
				}

				SearchToELIBJob job=new SearchToELIBJob(elibFile, processorTableModel);
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
			FileDialog dialog=new FileDialog(frame, "Select a MZML file", FileDialog.LOAD);
			dialog.setMultipleMode(true);
			dialog.setFilenameFilter(new SimpleFilenameFilter(".mzml", ".dia"));
			dialog.setVisible(true);
			if (dialog.getFiles()!=null) {
				for (File file : dialog.getFiles()) {
					getVisibleTab().getJob(file, processorTableModel);
				}
			}
		}
	}
}
