package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FileDialog;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.io.File;
import java.io.FilenameFilter;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFrame;
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
	
	public FileChooserPanel(File f, final String fileType, final FilenameFilter filter, boolean required) {
		super(new BorderLayout());
		this.required=required;
		fileLabel=new JLabel("Please select file...");
		chooseFile=new JButton("Edit");
		
		fileLabel.setBorder(BorderFactory.createLineBorder(Color.gray));
		
		chooseFile.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				JFrame frame = (JFrame)SwingUtilities.getRoot(FileChooserPanel.this);
				FileDialog dialog=new FileDialog(frame, "Select a "+fileType+" file", FileDialog.LOAD);
				if (file!=null) {
					dialog.setFile(file.getAbsolutePath());
				}
				dialog.setFilenameFilter(filter);
				dialog.setVisible(true);
				update(dialog.getFiles());
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
}
