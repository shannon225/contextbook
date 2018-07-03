package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dialog;
import java.awt.FileDialog;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.io.File;
import java.io.FilenameFilter;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

public class FileChooserPanel extends JPanel {
	private static final long serialVersionUID=1L;
	
	private File file=null;
	private final JLabel fileLabel;
	private final JButton chooseFile;
	private final JPanel top;
	private final boolean required;
	private final String fileType;
	private final FilenameFilter filter;
	private final boolean isLoad;

	public FileChooserPanel(File f, String fileTypeLocal, FilenameFilter filterLocal, boolean required) {
		this(f, fileTypeLocal, filterLocal, required, true);
	}
	public FileChooserPanel(File f, String fileTypeLocal, FilenameFilter filterLocal, boolean required, boolean isLoad) {
		super(new BorderLayout());
		this.setOpaque(true);
		this.setBackground(Color.white);
		
		this.isLoad=isLoad;
		this.fileType=fileTypeLocal;
		this.filter=filterLocal;
		this.required=required;
		fileLabel=new JLabel("Please select file...");
		chooseFile=new JButton("Edit");
		
		fileLabel.setBorder(BorderFactory.createLineBorder(Color.gray));
		
		chooseFile.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				askForFiles();
			}
		});
		
		top=new JPanel(new BorderLayout());
		top.setOpaque(true);
		
		this.add(top, BorderLayout.NORTH);
		
		top.add(new JLabel("<html><p style=\"font-size:10px; font-family: Helvetica, sans-serif\">"+fileType+": "), BorderLayout.WEST);
		top.add(fileLabel, BorderLayout.CENTER);
		top.add(chooseFile, BorderLayout.EAST);
		
		setToolTipText(fileLabel.getText());
		update(f);
	}

	public void askForFiles() {
		Component root=SwingUtilities.getRoot(FileChooserPanel.this);
		File[] files;
		if (root instanceof Frame) {
			files=getFiles(file, fileType, filter, (Frame)root, isLoad);
		} else {
			files=getFiles(file, fileType, filter, (Dialog)root, isLoad);
			
		}
		update(files);
	}
	
	@Override
	public String getToolTipText(MouseEvent event) {
		return fileLabel.getText();
	}
	
	public File getFile() {
		return file;
	}
	
	public void update(File... filename) {
		if (filename==null||filename.length==0||filename[0]==null) {
			file=null;
			fileLabel.setText("Please select file...");
			if (required) {
				top.setBackground(Color.pink);
			} else {
				top.setBackground(Color.lightGray);
			}
		} else {
			file=filename[0];
			fileLabel.setText(file.getName());
			top.setBackground(Color.white);
		}
	}

	public static File[] getFiles(File startingFile, String fileType, FilenameFilter filter, Dialog frame, boolean isLoad) {
		FileDialog dialog=new FileDialog(frame, "Select a "+fileType+" file", isLoad?FileDialog.LOAD:FileDialog.SAVE);
		if (startingFile!=null) {
			dialog.setFile(startingFile.getAbsolutePath());
		}
		dialog.setFilenameFilter(filter);
		dialog.setVisible(true);
		File[] files=dialog.getFiles();
		return files;
	}

	public static File[] getFiles(File startingFile, String fileType, FilenameFilter filter, Frame frame, boolean isLoad) {
		FileDialog dialog=new FileDialog(frame, "Select a "+fileType+" file", isLoad?FileDialog.LOAD:FileDialog.SAVE);
		if (startingFile!=null) {
			dialog.setFile(startingFile.getAbsolutePath());
		}
		dialog.setFilenameFilter(filter);
		dialog.setVisible(true);
		File[] files=dialog.getFiles();
		return files;
	}
}
