package edu.washington.gs.maccoss.encyclopedia.gui.framework;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.FileDialog;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.Optional;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.SwingUtilities;
import javax.swing.table.TableColumn;

import edu.washington.gs.maccoss.encyclopedia.gui.framework.library.EncyclopediaParametersPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.pecan.PecanParametersPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessorTableModel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LogConsole;
import edu.washington.gs.maccoss.encyclopedia.gui.general.MemoryMonitor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.ProgressRenderer;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class SearchPanel extends JPanel {
	private static final long serialVersionUID=1L;
	JobProcessorTableModel processorTableModel=new JobProcessorTableModel();
	
	private final JTabbedPane optionsTabs;
	
	public ParametersPanelInterface getVisibleTab() {
		return (ParametersPanelInterface)optionsTabs.getSelectedComponent();
	}
	
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
		JButton chooseFile=new JButton("Add MZML");
		chooseFile.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
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
		});
		
		JButton saveBlib=new JButton("Save BLIB");
		saveBlib.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
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
		});
		
		JButton saveElib=new JButton("Save Library");
		saveElib.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
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
        
		filesWrapper.add(table, BorderLayout.CENTER);
		
		split.setRightComponent(filesWrapper);
		split.setDividerLocation(600);

		this.add(split, BorderLayout.CENTER);
	}
}
