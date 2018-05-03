package edu.washington.gs.maccoss.encyclopedia.gui.framework;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Optional;
import java.util.StringTokenizer;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.SwingUtilities;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MSPReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.TraMLToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingWorkerProgress;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class SearchPanelUtilities {
	private static final ImageIcon convertDBIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/convertdb.png"));
	private static final ImageIcon fileAddIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/fileadd.png"));
	
	public static void combineELIBs(Component root) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Combine Libraries", true);
		
		final FileChooserPanel saveFileChooser=new FileChooserPanel(null, "Library", new SimpleFilenameFilter(".dlib", ".elib"), true, false);
		
		final JPanel choosers=new JPanel();
		choosers.setLayout(new BoxLayout(choosers, BoxLayout.Y_AXIS));
		choosers.add(new FileChooserPanel(null, "Add Library", new SimpleFilenameFilter(".dlib", ".elib"), false));
		
		JPanel organizer=new JPanel(new BorderLayout());
		organizer.add(choosers, BorderLayout.NORTH);
		JScrollPane scrollPane = new JScrollPane(organizer); 
		scrollPane.setPreferredSize(new Dimension(500, 400));
		JButton addChooserButton=new JButton("Add Additional Library Selector", fileAddIcon);
		addChooserButton.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				choosers.add(new FileChooserPanel(null, "Add Library", new SimpleFilenameFilter(".dlib", ".elib"), false), choosers.getComponentCount()-1);
				choosers.revalidate();
				choosers.repaint();
			}
		});

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(scrollPane);
		options.add(addChooserButton);
		options.add(saveFileChooser);
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				dialog.setVisible(false);
				dialog.dispose();
				
				final ArrayList<File> files=new ArrayList<>();
				for (Component c : choosers.getComponents()) {
					if (c instanceof FileChooserPanel) {
						File f=((FileChooserPanel)c).getFile();
						if (f!=null&&f.exists()) {
							files.add(f);
						}
					}
				}

				final File saveFile=saveFileChooser.getFile();
				
				if (files.size()>0&&saveFile!=null) {
					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame)SwingUtilities.getWindowAncestor(root), "Please wait...", "Reading Library Files") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							LibraryFile saveLibrary=new LibraryFile();
							saveLibrary.openFile();
							saveLibrary.dropIndices();
							
							for (File elibFile : files) {
								LibraryFile library=new LibraryFile();
								library.openFile(elibFile);
								final ArrayList<LibraryEntry> allEntries=library.getAllEntries(false,  new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()));
								saveLibrary.addEntries(allEntries);
								saveLibrary.addProteinsFromEntries(allEntries);
								Logger.logLine("Found "+allEntries.size()+" entries from "+elibFile.getName()+". Writing to ["+saveFile.getAbsolutePath()+"]...");
								library.close();
							}
							
							saveLibrary.createIndices();
							saveLibrary.saveAsFile(saveFile);
							
							saveLibrary.close();
							
							return Nothing.NOTHING;
						}
						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
					
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify at least one library file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
		dialog.setVisible(true);
	}
	
	public static void subsetELIB(Component root) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Subset Library", true);
		
		final FileChooserPanel elibFileChooser=new FileChooserPanel(null, "Library", new SimpleFilenameFilter(".dlib", ".elib"), true);
		final FileChooserPanel saveFileChooser=new FileChooserPanel(null, "Library", new SimpleFilenameFilter(".dlib", ".elib"), true, false);
		
		final JTextArea textArea = new JTextArea(25, 80);
		textArea.setFont(new Font("Monospaced", Font.PLAIN, 10));
		JScrollPane scrollPane = new JScrollPane(textArea); 
		textArea.setEditable(true);

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(elibFileChooser);
		options.add(saveFileChooser);
		options.add(new JLabel("Subset peptides:", JLabel.LEFT));
		options.add(scrollPane);
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				dialog.setVisible(false);
				dialog.dispose();

				final File elibFile=elibFileChooser.getFile();
				final File saveFile=saveFileChooser.getFile();
				final String text=textArea.getText();
				
				if (elibFile!=null&&elibFile.exists()&&saveFile!=null&&text!=null&&text.length()>0) {
					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame)SwingUtilities.getWindowAncestor(root), "Please wait...", "Reading Library File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							HashSet<String> targets=new HashSet<>();
							StringTokenizer st=new StringTokenizer(text);
							while (st.hasMoreTokens()) {
								targets.add(st.nextToken());
							}
							
							LibraryFile library=new LibraryFile();
							library.openFile(elibFile);
							

							LibraryFile saveLibrary=new LibraryFile();
							saveLibrary.openFile();
							
							ArrayList<LibraryEntry> toWrite=new ArrayList<>();
							for (LibraryEntry entry : library.getAllEntries(false, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()))) {
								if (targets.contains(entry.getPeptideSeq())) {
									toWrite.add(entry);
								}
							}
							Logger.logLine("Found "+toWrite.size()+" peptides from "+targets.size()+" target sequences. Writing to ["+saveFile.getAbsolutePath()+"]...");
							
							saveLibrary.dropIndices();
							saveLibrary.addEntries(toWrite);
							saveLibrary.addProteinsFromEntries(toWrite);
							saveLibrary.createIndices();
							saveLibrary.saveAsFile(saveFile);
							
							library.close();
							saveLibrary.close();
							
							return Nothing.NOTHING;
						}
						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
					
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify a library file and peptide sequences!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
		//dialog.setSize(500, 600);
		dialog.setVisible(true);
	}
	
	public static void convertBLIB(Component root, SearchParameters params) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
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

				final File blibFile=blibFileChooser.getFile();
				final File irtFile=iRTFileChooser.getFile();
				final File fastaFile=fastaFileChooser.getFile();
				
				if (blibFile!=null&&blibFile.exists()&&fastaFile!=null&&fastaFile.exists()) {
					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame)SwingUtilities.getWindowAncestor(root), "Please wait...", "Reading BLIB File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							BlibToLibraryConverter.convert(blibFile, Optional.ofNullable(irtFile), fastaFile, params);
							Logger.logLine("Finished reading "+blibFile.getName());
							return Nothing.NOTHING;
						}
						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
					
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
	
	public static void convertMSP(Component root, SearchParameters params) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Convert NIST MSP to Library", true);
		
		final FileChooserPanel mspFileChooser=new FileChooserPanel(null, "MSP", new SimpleFilenameFilter(".msp"), true);
		final FileChooserPanel fastaFileChooser=new FileChooserPanel(null, "FASTA", new SimpleFilenameFilter(".fas", ".fasta"), true);

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(mspFileChooser);
		options.add(fastaFileChooser);
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				dialog.setVisible(false);
				dialog.dispose();

				final File mspFile=mspFileChooser.getFile();
				final File fastaFile=fastaFileChooser.getFile();
				
				if (mspFile!=null&&mspFile.exists()&&fastaFile!=null&&fastaFile.exists()) {
					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame) SwingUtilities.getWindowAncestor(root), "Please wait...", "Reading MSP File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							MSPReader.convertMSP(mspFile, fastaFile);
							Logger.logLine("Finished reading "+mspFile.getName());
							return Nothing.NOTHING;
						}

						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
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
	
	public static void convertTRAML(Component root, SearchParameters params) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Convert TraML to Library", true);
		
		final FileChooserPanel tramlFileChooser=new FileChooserPanel(null, "TraML", new SimpleFilenameFilter(".traml"), true);
		final FileChooserPanel fastaFileChooser=new FileChooserPanel(null, "FASTA", new SimpleFilenameFilter(".fas", ".fasta"), true);

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(tramlFileChooser);
		options.add(fastaFileChooser);
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				dialog.setVisible(false);
				dialog.dispose();

				final File tramlFile=tramlFileChooser.getFile();
				final File fastaFile=fastaFileChooser.getFile();
				
				if (tramlFile!=null&&tramlFile.exists()&&fastaFile!=null&&fastaFile.exists()) {
					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame) SwingUtilities.getWindowAncestor(root), "Please wait...", "Reading TraML File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							TraMLToLibraryConverter.convertTraML(tramlFile, fastaFile, params.getAAConstants());
							Logger.logLine("Finished reading "+tramlFile.getName());
							return Nothing.NOTHING;
						}

						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify a TraML and a FASTA file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
}
