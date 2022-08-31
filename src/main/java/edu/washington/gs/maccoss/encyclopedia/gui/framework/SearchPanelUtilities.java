package edu.washington.gs.maccoss.encyclopedia.gui.framework;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Desktop;
import java.awt.FileDialog;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Optional;
import java.util.StringTokenizer;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTextArea;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaToPrositCSVParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryToBlibConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MS2PIPReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MSPReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MaxquantMSMSConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.OpenSwathTSVToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SpectronautCSVToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.TraMLSAXToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filewriters.LibraryUtilities;
import edu.washington.gs.maccoss.encyclopedia.filewriters.MS2PIPWriter;
import edu.washington.gs.maccoss.encyclopedia.filewriters.MSPWriter;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PrositCSVWriter;
import edu.washington.gs.maccoss.encyclopedia.filewriters.StripeFileMerger;
import edu.washington.gs.maccoss.encyclopedia.filewriters.StripeFileTrimmer;
import edu.washington.gs.maccoss.encyclopedia.gui.general.AboutDialog;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserList;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LabeledComponent;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingWorkerProgress;
import edu.washington.gs.maccoss.encyclopedia.jobs.WorkerJob;
import edu.washington.gs.maccoss.encyclopedia.jobs.XMLDriverFactory;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.StringUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLObject;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class SearchPanelUtilities {
	private static final ImageIcon convertDBIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/convertdb.png"));
	private static final ImageIcon fileAddIcon=new ImageIcon(SearchPanel.class.getClassLoader().getResource("images/fileadd.png"));

	public static void preprocessMZMLs(Component root, SearchParameters params) {

		FileDialog dialog=new FileDialog((JFrame)null, "mzML files", FileDialog.LOAD);
		dialog.setMultipleMode(true);
		dialog.setFilenameFilter(new SimpleFilenameFilter(".mzML", StripeFile.DIA_EXTENSION));
		dialog.setVisible(true);
		File[] featureFiles=dialog.getFiles();

		if (featureFiles!=null) {
			if (featureFiles.length==1&&featureFiles[0].isDirectory()) {
				featureFiles=featureFiles[0].listFiles();
			}
			
			final File[] filesToLoad=featureFiles;
			
			SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame)SwingUtilities.getWindowAncestor(root), "Please wait...", "Reading mzML Files") {
				/**
				 * FIXME ADD PROGRESS
				 * @param root
				 * @param params
				 */
				@Override
				protected Nothing doInBackgroundForReal() throws Exception {
					for (int i=0; i<filesToLoad.length; i++) {
						System.out.println("Processing "+filesToLoad[i]);
						StripeFileInterface file=StripeFileGenerator.getFile(filesToLoad[i], params);
						if (filesToLoad[i].getName().endsWith(StripeFile.DIA_EXTENSION) && file instanceof StripeFile) {
							((StripeFile)file).saveFile();
						}
					}
					return Nothing.NOTHING;
				}
				@Override
				protected void doneForReal(Nothing t) {
				}
			};
			worker.execute();
		}
	}
	
	public static void convertELIBtoOpenSWATH(Component root, SearchParameters params) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Convert Library to OpenSWATH", true);

		final FileChooserPanel elibFileChooser=new FileChooserPanel(null, "Library", new SimpleFilenameFilter(".dlib", ".elib"), true, true);

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(elibFileChooser);
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				final File elibFile=elibFileChooser.getFile();
				String absolutePath=elibFile.getAbsolutePath();
				File tsvFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+".tsv");

				if (elibFile!=null&&elibFile.exists()) {
					dialog.setVisible(false);
					dialog.dispose();

					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame)SwingUtilities.getWindowAncestor(root), "Please wait...", "Reading Library File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							OpenSwathTSVToLibraryConverter.convertToOpenSwathTSV(params, elibFile, tsvFile);
							return Nothing.NOTHING;
						}
						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify an ELIB or DLIB library file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
	
	public static void convertELIBtoMSP(Component root, SearchParameters params) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Convert Library to NIST MSP", true);

		final FileChooserPanel elibFileChooser=new FileChooserPanel(null, "Library", new SimpleFilenameFilter(".dlib", ".elib"), true, true);
		final JCheckBox filterToKnownIonsBox=new JCheckBox("Filter out unknown ions");
		filterToKnownIonsBox.setSelected(false);

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(elibFileChooser);
		options.add(filterToKnownIonsBox);
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				final File elibFile=elibFileChooser.getFile();
				String absolutePath=elibFile.getAbsolutePath();
				File mspFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+".msp");
				final boolean filterToKnownIons=filterToKnownIonsBox.isSelected();

				if (elibFile!=null&&elibFile.exists()) {
					dialog.setVisible(false);
					dialog.dispose();

					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame)SwingUtilities.getWindowAncestor(root), "Please wait...", "Reading Library File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							LibraryFile library=new LibraryFile();
							library.openFile(elibFile);
							MSPWriter.writeMSP(mspFile, library, filterToKnownIons, params);
							library.close();
							return Nothing.NOTHING;
						}
						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify an ELIB or DLIB library file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
	
	public static void convertELIBtoBLIB(Component root) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Convert Library to BLIB", true);

		final FileChooserPanel elibFileChooser=new FileChooserPanel(null, "Library", new SimpleFilenameFilter(".dlib", ".elib"), true, true);

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(elibFileChooser);
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				final File elibFile=elibFileChooser.getFile();
				String absolutePath=elibFile.getAbsolutePath();
				File blibFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+BlibFile.BLIB);
				
				if (elibFile!=null&&elibFile.exists()) {
					dialog.setVisible(false);
					dialog.dispose();

					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame)SwingUtilities.getWindowAncestor(root), "Please wait...", "Reading Library File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							LibraryToBlibConverter.convert(elibFile, blibFile);
							return Nothing.NOTHING;
						}
						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
					
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify an ELIB or DLIB library file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
	
	public static void saveDriverFile(Component root, final SearchParameters params, final JobProcessor processor) {
		JFrame frame = (JFrame)SwingUtilities.getRoot(root);

		ArrayList<SearchJobData> jobData=new ArrayList<SearchJobData>();
		for (WorkerJob job : processor.getQueue()) {
			if (job instanceof SearchJob) {
				jobData.add(((SearchJob)job).getSearchData());
			}
		}
		if (jobData.size()<1) {
			JOptionPane.showMessageDialog(frame, "Please queue at least one RAW file first!");
			
		} else {
			FileDialog dialog=new FileDialog(frame, "Select an XML file", FileDialog.SAVE);
			SimpleFilenameFilter filter = new SimpleFilenameFilter(XMLDriverFactory.DRIVER_XML_EXTENSION);
			dialog.setFilenameFilter(filter);
			dialog.setVisible(true);
			File[] files=dialog.getFiles();
			if (files!=null&&files.length>0) {
				ArrayList<WorkerJob> queue=processor.getQueue();
				File file=files[0];
				if (!filter.accept(file.getName())) {
					file=new File(file.getParent(), file.getName()+XMLDriverFactory.DRIVER_XML_EXTENSION);
				}

				XMLDriverFactory.writeXML(queue, file);
			}
		}
	}
	
	public static void extractSampleSpecificDLIBs(Component root, final SearchParameters params, final JobProcessor processor) {
		JFrame frame = (JFrame)SwingUtilities.getRoot(root);

		ArrayList<SearchJobData> jobData=new ArrayList<SearchJobData>();
		for (WorkerJob job : processor.getQueue()) {
			if (job instanceof SearchJob) {
				jobData.add(((SearchJob)job).getSearchData());
			}
		}
		if (jobData.size()<2) {
			JOptionPane.showMessageDialog(frame, "Please queue at least two RAW files first!");
			
		} else {
			final JDialog dialog=new JDialog(frame, "Extact Sample-Specific Libraries from ELIB", true);
	
			final FileChooserPanel diaFileChooser=new FileChooserPanel(null, "Single-injection DIA example", new SimpleFilenameFilter(MzmlToDIAConverter.MZML_EXTENSION, StripeFile.DIA_EXTENSION), true, true);
			final FileChooserPanel saveDirFileChooser=new FileChooserPanel(null, "Save Directory", new FilenameFilter() {
				@Override
				public boolean accept(File dir, String name) {
					if (dir.exists()&&!dir.isDirectory()) {
						return false;
					}
					return true;
				}
			}, true, false);
	
			JPanel options=new JPanel();
			options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
			options.add(diaFileChooser);
			options.add(saveDirFileChooser);
			
			JPanel buttons=new JPanel();
			buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
			JButton okButton=new JButton("OK");
			okButton.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					final File saveDir=saveDirFileChooser.getFile();
					final File exampleDIAFile=diaFileChooser.getFile(); // can be missing (null)
					if (saveDir!=null) {
						dialog.setVisible(false);
						dialog.dispose();
	
						Logger.logLine("Added extact sample-specific libraries into ["+saveDir.getAbsolutePath()+"]");
						processor.addJob(new CombineELIBsAndExtractGroupSpecificLibrariesJob(saveDir, Optional.of(exampleDIAFile), params, processor));
					} else {
						JOptionPane.showMessageDialog(frame, "You must specify an ELIB or DLIB library file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
	}
	
	public static void combineELIBs(Component root, final JobProcessor processor) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Combine Libraries", true);
		
		final FileChooserPanel saveFileChooser=new FileChooserPanel(null, "Library", new SimpleFilenameFilter(".dlib"), true, false);

		final FileChooserList choosers=new FileChooserList("Library File (.dlib or .elib)", new SimpleFilenameFilter(".dlib", ".elib"));
		
		final JCheckBox rtAlignBox=new JCheckBox("RT align samples");
		final JCheckBox removeDuplicatesBox=new JCheckBox("Remove duplicates");
		final JCheckBox higherScoresAreBetterBox=new JCheckBox("Higher scores are better");
		rtAlignBox.setSelected(false);
		removeDuplicatesBox.setSelected(true);
		higherScoresAreBetterBox.setSelected(false);

		removeDuplicatesBox.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
				higherScoresAreBetterBox.setEnabled(removeDuplicatesBox.isSelected());
			}
		});

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(choosers);
		options.add(saveFileChooser);
		options.add(rtAlignBox);
		options.add(removeDuplicatesBox);
		options.add(higherScoresAreBetterBox);
		
		choosers.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				ArrayList<File> filelist=choosers.getFiles();
				if (filelist.size()>0) {
					// choose the last file's directory to place the new file
					File directory=filelist.get(filelist.size()-1).getParentFile();
					
					String[] names=new String[filelist.size()];
					for (int i = 0; i < names.length; i++) {
						names[i]=filelist.get(i).getName();
					}
					String filename=StringUtils.getCommonName(names, "combined");
					if (filename.toLowerCase().endsWith(LibraryFile.ELIB)) {
						filename=filename.substring(0, filename.length()-LibraryFile.ELIB.length())+LibraryFile.DLIB;
					} else if (!filename.toLowerCase().endsWith(LibraryFile.DLIB)) {
						filename=filename+LibraryFile.DLIB;
					}
					
					File savefile=new File(directory, filename);
					saveFileChooser.update(savefile);
				}
			}
		});
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				final ArrayList<File> files=new ArrayList<>();
				for (File f : choosers.getFiles()) {
					if (f!=null&&f.exists()) {
						files.add(f);
					}
				}

				final File saveFile=saveFileChooser.getFile();
				final boolean higherScoresAreBetter=higherScoresAreBetterBox.isSelected();
				final boolean rtAlign=rtAlignBox.isSelected();
				final boolean removeDuplicates=removeDuplicatesBox.isSelected();
				
				if (files.size()>0&&saveFile!=null) {
					dialog.setVisible(false);
					dialog.dispose();

					WorkerJob job=new WorkerJob() {
						@Override
						public String getJobTitle() {
							return "Merge raw files into "+saveFile.getName();
						}

						@Override
						public void runJob(ProgressIndicator progress) throws Exception {
							LibraryUtilities.mergeLibraries(progress, files, saveFile, rtAlign, removeDuplicates, higherScoresAreBetter);
						}
					};
					processor.addJob(job);
					
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
	
	public static void combineMZMLs(Component root, final JobProcessor processor, SearchParameters parameters) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Combine DIA or mzML Gas Phase Fractions", true);
		
		final FileChooserPanel saveFileChooser=new FileChooserPanel(null, "DIA File", new SimpleFilenameFilter(".dia"), true, false);
		
		final FileChooserList choosers=new FileChooserList("GPF DIA/mzML File", new SimpleFilenameFilter(".mzML", ".dia"));
		

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(choosers);
		options.add(saveFileChooser);
		
		choosers.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				ArrayList<File> filelist=choosers.getFiles();
				if (filelist.size()>0) {
					// choose the last file's directory to place the new file
					File directory=filelist.get(filelist.size()-1).getParentFile();
					
					String[] names=new String[filelist.size()];
					for (int i = 0; i < names.length; i++) {
						names[i]=filelist.get(i).getName();
					}
					String filename=StringUtils.getCommonName(names, "combined");
					if (filename.toLowerCase().endsWith(".mzml")) {
						filename=filename.substring(0, filename.length()-5)+StripeFile.DIA_EXTENSION;
					} else if (!filename.toLowerCase().endsWith(StripeFile.DIA_EXTENSION)) {
						filename=filename+StripeFile.DIA_EXTENSION;
					}
					
					File savefile=new File(directory, filename);
					saveFileChooser.update(savefile);
				}
			}
		});
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				final ArrayList<File> files=new ArrayList<>();
				for (File f : choosers.getFiles()) {
					if (f!=null&&f.exists()) {
						files.add(f);
					}
				}

				final File saveFile=saveFileChooser.getFile();
				
				if (files.size()>0&&saveFile!=null) {
					dialog.setVisible(false);
					dialog.dispose();
					
					WorkerJob job=new WorkerJob() {
						@Override
						public String getJobTitle() {
							return "Merge raw files into "+saveFile.getName();
						}

						@Override
						public void runJob(ProgressIndicator progress) throws Exception {
							StripeFileMerger.merge(progress, files.toArray(new File[files.size()]), saveFile, parameters);
						}
					};
					processor.addJob(job);
					
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify at least one DIA/mzML file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
	
	public static void subsetDIA(Component root, SearchParameters params) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Subset Raw File", true);
		
		final FileChooserPanel diaFileChooser=new FileChooserPanel(null, "Starting Raw File (mzML, DIA)", new SimpleFilenameFilter(".dia", ".mzML"), true);
		final FileChooserPanel saveFileChooser=new FileChooserPanel(null, "New Raw File (DIA)", new SimpleFilenameFilter(".dia", ".mzML"), true, false);
		
		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(diaFileChooser);
		options.add(saveFileChooser);
		
		final SpinnerModel minRT=new SpinnerNumberModel(0.0, 0.0, 1000.0, 1.0);
		final SpinnerModel maxRT=new SpinnerNumberModel(1000.0, 0.0, 1000.0, 1.0);
		JPanel rtRange=new JPanel(new FlowLayout());
		rtRange.setOpaque(true);
		rtRange.setBackground(Color.white);
		rtRange.add(new JSpinner(minRT));
		rtRange.add(new JLabel("<html><p style=\"font-size:10px; font-family: Helvetica, sans-serif\"> to "));
		rtRange.add(new JSpinner(maxRT));

		final SpinnerModel minMZ=new SpinnerNumberModel(0.0, 0.0, 1000.0, 1.0);
		final SpinnerModel maxMZ=new SpinnerNumberModel(10000.0, 0.0, 10000.0, 1.0);
		JPanel mzRange=new JPanel(new FlowLayout());
		mzRange.setOpaque(true);
		mzRange.setBackground(Color.white);
		mzRange.add(new JSpinner(minMZ));
		mzRange.add(new JLabel("<html><p style=\"font-size:10px; font-family: Helvetica, sans-serif\"> to "));
		mzRange.add(new JSpinner(maxMZ));

		options.add(new LabeledComponent("Precursor Range (m/z)", mzRange));
		options.add(new LabeledComponent("Retention Time Range (min)", rtRange));
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				final File diaFile=diaFileChooser.getFile();
				final File saveFile=saveFileChooser.getFile();

				final double mzMin=((Number)minMZ.getValue()).doubleValue();
				final double mzMax=((Number)maxMZ.getValue()).doubleValue();
				if (mzMax<mzMin) {
					JOptionPane.showMessageDialog(frame, "Precursor m/z maximum must be higher than the minimum!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
					return;
				}

				final float rtMinSec=60f*((Number)minRT.getValue()).floatValue();
				final float rtMaxSec=60f*((Number)maxRT.getValue()).floatValue();
				if (rtMaxSec<rtMinSec) {
					JOptionPane.showMessageDialog(frame, "Retention time maximum must be higher than the minimum!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
					return;
				}
				
				if (diaFile!=null&&diaFile.exists()&&saveFile!=null) {
					dialog.setVisible(false);
					dialog.dispose();

					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame)SwingUtilities.getWindowAncestor(root), "Please wait...", "Reading Library File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							StripeFileTrimmer.trim(diaFile, saveFile, new Range(mzMin, mzMax), new Range(rtMinSec, rtMaxSec), params);
							
							return Nothing.NOTHING;
						}
						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
					
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify a raw file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
	
	public static void subsetELIB(Component root) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Subset Library", true);
		
		final FileChooserPanel elibFileChooser=new FileChooserPanel(null, "Starting Library", new SimpleFilenameFilter(".dlib", ".elib"), true);
		final FileChooserPanel saveFileChooser=new FileChooserPanel(null, "New Library", new SimpleFilenameFilter(".dlib", ".elib"), true, false);
		
		final JTextArea textArea = new JTextArea(25, 80);
		textArea.setFont(new Font("Monospaced", Font.PLAIN, 10));
		JScrollPane scrollPane = new JScrollPane(textArea); 
		textArea.setEditable(true);

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(elibFileChooser);
		options.add(saveFileChooser);
		
		final SpinnerModel minRT=new SpinnerNumberModel(0.0, 0.0, 1000.0, 1.0);
		final SpinnerModel maxRT=new SpinnerNumberModel(1000.0, 0.0, 1000.0, 1.0);
		JPanel rtRange=new JPanel(new FlowLayout());
		rtRange.setOpaque(true);
		rtRange.setBackground(Color.white);
		rtRange.add(new JSpinner(minRT));
		rtRange.add(new JLabel("<html><p style=\"font-size:10px; font-family: Helvetica, sans-serif\"> to "));
		rtRange.add(new JSpinner(maxRT));

		final SpinnerModel minMZ=new SpinnerNumberModel(0.0, 0.0, 1000.0, 1.0);
		final SpinnerModel maxMZ=new SpinnerNumberModel(10000.0, 0.0, 10000.0, 1.0);
		JPanel mzRange=new JPanel(new FlowLayout());
		mzRange.setOpaque(true);
		mzRange.setBackground(Color.white);
		mzRange.add(new JSpinner(minMZ));
		mzRange.add(new JLabel("<html><p style=\"font-size:10px; font-family: Helvetica, sans-serif\"> to "));
		mzRange.add(new JSpinner(maxMZ));

		options.add(new LabeledComponent("Precursor Range (m/z)", mzRange));
		options.add(new LabeledComponent("Rention Time Range (min)", rtRange));
		
		options.add(new JLabel("Subset peptides/accessions (requires exact matches):", JLabel.LEFT));
		options.add(scrollPane);
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				final File elibFile=elibFileChooser.getFile();
				final File saveFile=saveFileChooser.getFile();

				final double mzMin=((Number)minMZ.getValue()).doubleValue();
				final double mzMax=((Number)maxMZ.getValue()).doubleValue();
				if (mzMax<mzMin) {
					JOptionPane.showMessageDialog(frame, "Precursor m/z maximum must be higher than the minimum!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
					return;
				}

				final float rtMinSec=60f*((Number)minRT.getValue()).floatValue();
				final float rtMaxSec=60f*((Number)maxRT.getValue()).floatValue();
				if (rtMaxSec<rtMinSec) {
					JOptionPane.showMessageDialog(frame, "Retention time maximum must be higher than the minimum!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
					return;
				}
				
				final String text=textArea.getText();

				final HashSet<String> targets=new HashSet<>();
				if (text!=null&&text.length()>0) {
					StringTokenizer st=new StringTokenizer(text);
					while (st.hasMoreTokens()) {
						targets.add(st.nextToken());
					}
				}
				
				if (elibFile!=null&&elibFile.exists()&&saveFile!=null) {
					dialog.setVisible(false);
					dialog.dispose();

					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame)SwingUtilities.getWindowAncestor(root), "Please wait...", "Reading Library File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							
							LibraryFile library=new LibraryFile();
							library.openFile(elibFile);
							
							LibraryUtilities.subsetLibrary(saveFile, rtMinSec, rtMaxSec, mzMin, mzMax, targets, library);
							
							library.close();
							
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
		//final JCheckBox higherScoresAreBetterBox=new JCheckBox("Higher scores are better");
		//higherScoresAreBetterBox.setSelected(false);

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(blibFileChooser);
		options.add(iRTFileChooser);
		options.add(fastaFileChooser);
		//options.add(higherScoresAreBetterBox);
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				final File blibFile=blibFileChooser.getFile();
				final File irtFile=iRTFileChooser.getFile();
				final File fastaFile=fastaFileChooser.getFile();
				//final boolean higherScoresAreBetter=higherScoresAreBetterBox.isSelected();
				
				if (blibFile!=null&&blibFile.exists()&&fastaFile!=null&&fastaFile.exists()) {
					dialog.setVisible(false);
					dialog.dispose();

					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame)SwingUtilities.getWindowAncestor(root), "Please wait...", "Reading BLIB File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							BlibToLibraryConverter.convert(blibFile, Optional.ofNullable(irtFile), fastaFile, false, params);
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
		dialog.setSize(500, 250);
		dialog.setVisible(true);
	}
	
	public static void convertFastaForProsit(Component root) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Convert FASTA to Prosit CSV", true);
		
		
		final FileChooserPanel fastaFileChooser=new FileChooserPanel(null, "FASTA", new SimpleFilenameFilter(".fas", ".fasta"), true);

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(fastaFileChooser);

		final SpinnerModel defaultNCESpinner=new SpinnerNumberModel(FastaToPrositCSVParameters.DEFAULT_DEFAULT_NCE, FastaToPrositCSVParameters.MIN_DEFAULT_NCE, FastaToPrositCSVParameters.MAX_DEFAULT_NCE, 1);
		final SpinnerModel defaultChargeSpinner=new SpinnerNumberModel(FastaToPrositCSVParameters.DEFAULT_DEFAULT_CHARGE, FastaToPrositCSVParameters.MIN_DEFAULT_CHARGE, FastaToPrositCSVParameters.MAX_DEFAULT_CHARGE, 1);
		final SpinnerModel minChargeSpinner=new SpinnerNumberModel(FastaToPrositCSVParameters.DEFAULT_MIN_CHARGE, FastaToPrositCSVParameters.MIN_CHARGE, FastaToPrositCSVParameters.MAX_CHARGE, 1);
		final SpinnerModel maxChargeSpinner=new SpinnerNumberModel(FastaToPrositCSVParameters.DEFAULT_MAX_CHARGE, FastaToPrositCSVParameters.MIN_CHARGE, FastaToPrositCSVParameters.MAX_CHARGE, 1);
		final SpinnerModel maxMissedCleavageSpinner=new SpinnerNumberModel(FastaToPrositCSVParameters.DEFAULT_MAX_MISSED_CLEAVAGE, FastaToPrositCSVParameters.MIN_MAX_MISSED_CLEAVAGE, FastaToPrositCSVParameters.MAX_MAX_MISSED_CLEAVAGE, 1);
		final SpinnerModel minMzSpinner=new SpinnerNumberModel(FastaToPrositCSVParameters.DEFAULT_MIN_MZ, FastaToPrositCSVParameters.MIN_MZ, FastaToPrositCSVParameters.MAX_MZ, 0.1);
		final SpinnerModel maxMzSpinner=new SpinnerNumberModel(FastaToPrositCSVParameters.DEFAULT_MAX_MZ, FastaToPrositCSVParameters.MIN_MZ, FastaToPrositCSVParameters.MAX_MZ, 0.1);
		final JComboBox<String> enzymeBox=new JComboBox<String>(new String[] {"Trypsin", "Glu-C", "Lys-C", "Arg-C", "Asp-N", "Lys-N", "CNBr", "Chymotrypsin", "Pepsin A", "No Enzyme"});
		
		JPanel chargeRange=new JPanel(new FlowLayout());
		chargeRange.setOpaque(true);
		chargeRange.setBackground(Color.white);
		chargeRange.add(new JSpinner(minChargeSpinner));
		chargeRange.add(new JLabel("<html><p style=\"font-size:10px; font-family: Helvetica, sans-serif\"> to "));
		chargeRange.add(new JSpinner(maxChargeSpinner));
		//options.add(new LabeledComponent("Enzyme", enzymeBox)); // FIXME add prosit enzymes
		options.add(new LabeledComponent("Charge range", chargeRange));
		options.add(new LabeledComponent("Maximum Missed Cleavage", new JSpinner(maxMissedCleavageSpinner)));

		JPanel mzRange=new JPanel(new FlowLayout());
		mzRange.setOpaque(true);
		mzRange.setBackground(Color.white);
		mzRange.add(new JSpinner(minMzSpinner));
		mzRange.add(new JLabel("<html><p style=\"font-size:10px; font-family: Helvetica, sans-serif\"> to "));
		mzRange.add(new JSpinner(maxMzSpinner));
		options.add(new LabeledComponent("m/z range", mzRange));

		options.add(new LabeledComponent("Default NCE", "NCE is the Normalized Collision Energy for Thermo Fusion-class instruments. If you use QEs, add 6 to your NCE. If you use ToFs, use NCE=33.", new JSpinner(defaultNCESpinner)));
		options.add(new LabeledComponent("Default Charge", new JSpinner(defaultChargeSpinner)));
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				final File fastaFile=fastaFileChooser.getFile();
				byte defaultNCE=((Number)defaultNCESpinner.getValue()).byteValue();
				byte defaultCharge=((Number)defaultChargeSpinner.getValue()).byteValue();
				byte minCharge=((Number)minChargeSpinner.getValue()).byteValue();
				byte maxCharge=((Number)maxChargeSpinner.getValue()).byteValue();
				int maxMissedCleavages=((Number)maxMissedCleavageSpinner.getValue()).byteValue();
				double minimumMz=((Number)minMzSpinner.getValue()).doubleValue();
				double maximumMz=((Number)maxMzSpinner.getValue()).doubleValue();
				DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme((String)enzymeBox.getSelectedItem());
				
				if (fastaFile!=null&&fastaFile.exists()) {
					dialog.setVisible(false);
					dialog.dispose();
					
					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame) SwingUtilities.getWindowAncestor(root), "Please wait...", "Creating Prosit CSV File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							PrositCSVWriter.writeCSV(fastaFile, enzyme, defaultNCE, defaultCharge, minCharge, maxCharge, maxMissedCleavages, new Range(minimumMz, maximumMz), false);
							return Nothing.NOTHING;
						}

						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify a FASTA file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
		final JLabel text=new JLabel(AboutDialog.citationIcon);
		text.setText("<html><p style=\"font-size:10px; font-family: Helvetica, sans-serif\">"+"This function will <i>in silico</i> digest peptides from your FASTA and create an input file for Prosit. If you use this feature, please cite <a href=\"https://www.nature.com/articles/s41467-020-15346-1\">Searle et al, 2020</a>.");
		text.setBackground(Color.WHITE);
		text.setOpaque(true);
		text.addMouseListener(new MouseListener() {
			public void mouseClicked(MouseEvent e) {
			    try {
			         
			        Desktop.getDesktop().browse(new URI("https://www.nature.com/articles/s41467-020-15346-1"));
			         
			    } catch (IOException | URISyntaxException e1) {
			        e1.printStackTrace();
			    }
			}
			@Override
			public void mouseReleased(MouseEvent e) {
			}
			
			@Override
			public void mousePressed(MouseEvent e) {
			}
			
			@Override
			public void mouseExited(MouseEvent e) {
			}
			
			@Override
			public void mouseEntered(MouseEvent e) {
			}
		});

		mainpane.add(text, BorderLayout.NORTH);
		mainpane.add(options, BorderLayout.CENTER);
		mainpane.add(buttons, BorderLayout.SOUTH);
		mainpane.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10), BorderFactory.createTitledBorder("Parameters:")));
		
		dialog.getContentPane().add(mainpane, BorderLayout.CENTER);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		dialog.pack(); 
		dialog.setSize(500, 350);
		dialog.setVisible(true);
	}
	
	public static void convertFastaForMS2PIP(Component root) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Convert FASTA to MS2PIP PEPREC", true);
		
		
		final FileChooserPanel fastaFileChooser=new FileChooserPanel(null, "FASTA", new SimpleFilenameFilter(".fas", ".fasta"), true);

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(fastaFileChooser);

		final SpinnerModel minChargeSpinner=new SpinnerNumberModel(FastaToPrositCSVParameters.DEFAULT_MIN_CHARGE, FastaToPrositCSVParameters.MIN_CHARGE, FastaToPrositCSVParameters.MAX_CHARGE, 1);
		final SpinnerModel maxChargeSpinner=new SpinnerNumberModel(FastaToPrositCSVParameters.DEFAULT_MAX_CHARGE, FastaToPrositCSVParameters.MIN_CHARGE, FastaToPrositCSVParameters.MAX_CHARGE, 1);
		final SpinnerModel maxMissedCleavageSpinner=new SpinnerNumberModel(FastaToPrositCSVParameters.DEFAULT_MAX_MISSED_CLEAVAGE, FastaToPrositCSVParameters.MIN_MAX_MISSED_CLEAVAGE, FastaToPrositCSVParameters.MAX_MAX_MISSED_CLEAVAGE, 1);
		final SpinnerModel minMzSpinner=new SpinnerNumberModel(FastaToPrositCSVParameters.DEFAULT_MIN_MZ, FastaToPrositCSVParameters.MIN_MZ, FastaToPrositCSVParameters.MAX_MZ, 0.1);
		final SpinnerModel maxMzSpinner=new SpinnerNumberModel(FastaToPrositCSVParameters.DEFAULT_MAX_MZ, FastaToPrositCSVParameters.MIN_MZ, FastaToPrositCSVParameters.MAX_MZ, 0.1);
		final JComboBox<String> enzymeBox=new JComboBox<String>(new String[] {"Trypsin", "Glu-C", "Lys-C", "Arg-C", "Asp-N", "Lys-N", "CNBr", "Chymotrypsin", "Pepsin A", "No Enzyme"});
		
		JPanel chargeRange=new JPanel(new FlowLayout());
		chargeRange.setOpaque(true);
		chargeRange.setBackground(Color.white);
		chargeRange.add(new JSpinner(minChargeSpinner));
		chargeRange.add(new JLabel("<html><p style=\"font-size:10px; font-family: Helvetica, sans-serif\"> to "));
		chargeRange.add(new JSpinner(maxChargeSpinner));
		options.add(new LabeledComponent("Enzyme", enzymeBox)); // FIXME add prosit enzymes
		options.add(new LabeledComponent("Charge range", chargeRange));
		options.add(new LabeledComponent("Maximum Missed Cleavage", new JSpinner(maxMissedCleavageSpinner)));

		JPanel mzRange=new JPanel(new FlowLayout());
		mzRange.setOpaque(true);
		mzRange.setBackground(Color.white);
		mzRange.add(new JSpinner(minMzSpinner));
		mzRange.add(new JLabel("<html><p style=\"font-size:10px; font-family: Helvetica, sans-serif\"> to "));
		mzRange.add(new JSpinner(maxMzSpinner));
		options.add(new LabeledComponent("m/z range", mzRange));

		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				final File fastaFile=fastaFileChooser.getFile();
				byte minCharge=((Number)minChargeSpinner.getValue()).byteValue();
				byte maxCharge=((Number)maxChargeSpinner.getValue()).byteValue();
				int maxMissedCleavages=((Number)maxMissedCleavageSpinner.getValue()).byteValue();
				double minimumMz=((Number)minMzSpinner.getValue()).doubleValue();
				double maximumMz=((Number)maxMzSpinner.getValue()).doubleValue();
				DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme((String)enzymeBox.getSelectedItem());
				
				if (fastaFile!=null&&fastaFile.exists()) {
					dialog.setVisible(false);
					dialog.dispose();
					
					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame) SwingUtilities.getWindowAncestor(root), "Please wait...", "Creating MS2PIP PEPREC File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							MS2PIPWriter.writeMS2PIP(fastaFile, enzyme, minCharge, maxCharge, maxMissedCleavages, new Range(minimumMz, maximumMz), false);
							return Nothing.NOTHING;
						}

						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify a FASTA file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
		final JLabel text=new JLabel(AboutDialog.citationIcon);
		text.setText("<html><p style=\"font-size:10px; font-family: Helvetica, sans-serif\">"+"This function will <i>in silico</i> digest peptides from your FASTA and create an input file for MS2PIP. If you use this feature, please cite <a href=\"https://www.nature.com/articles/s41467-020-15346-1\">Searle et al, 2020</a>.");
		text.setBackground(Color.WHITE);
		text.setOpaque(true);
		text.addMouseListener(new MouseListener() {
			public void mouseClicked(MouseEvent e) {
			    try {
			         
			        Desktop.getDesktop().browse(new URI("https://www.nature.com/articles/s41467-020-15346-1"));
			         
			    } catch (IOException | URISyntaxException e1) {
			        e1.printStackTrace();
			    }
			}
			@Override
			public void mouseReleased(MouseEvent e) {
			}
			
			@Override
			public void mousePressed(MouseEvent e) {
			}
			
			@Override
			public void mouseExited(MouseEvent e) {
			}
			
			@Override
			public void mouseEntered(MouseEvent e) {
			}
		});

		mainpane.add(text, BorderLayout.NORTH);
		mainpane.add(options, BorderLayout.CENTER);
		mainpane.add(buttons, BorderLayout.SOUTH);
		mainpane.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10), BorderFactory.createTitledBorder("Parameters:")));
		
		dialog.getContentPane().add(mainpane, BorderLayout.CENTER);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		dialog.pack(); 
		dialog.setSize(500, 350);
		dialog.setVisible(true);
	}
	
	public static void convertLibraryForProsit(Component root) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Convert Library to Prosit CSV", true);
		
		
		final FileChooserPanel libraryFileChooser=new FileChooserPanel(null, "Library", new SimpleFilenameFilter(".dlib", ".elib"), true);

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(libraryFileChooser);

		final SpinnerModel defaultNCESpinner=new SpinnerNumberModel(FastaToPrositCSVParameters.DEFAULT_DEFAULT_NCE, FastaToPrositCSVParameters.MIN_DEFAULT_NCE, FastaToPrositCSVParameters.MAX_DEFAULT_NCE, 1);
		final SpinnerModel defaultChargeSpinner=new SpinnerNumberModel(FastaToPrositCSVParameters.DEFAULT_DEFAULT_CHARGE, FastaToPrositCSVParameters.MIN_DEFAULT_CHARGE, FastaToPrositCSVParameters.MAX_DEFAULT_CHARGE, 1);
		
		options.add(new LabeledComponent("Default NCE", "NCE is the Normalized Collision Energy for Thermo Fusion-class instruments. If you use QEs, add 6 to your NCE. If you use ToFs, use NCE=33.", new JSpinner(defaultNCESpinner)));
		options.add(new LabeledComponent("Default Charge", new JSpinner(defaultChargeSpinner)));
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				final File libraryFile=libraryFileChooser.getFile();
				byte defaultNCE=((Number)defaultNCESpinner.getValue()).byteValue();
				byte defaultCharge=((Number)defaultChargeSpinner.getValue()).byteValue();
				
				if (libraryFile!=null&&libraryFile.exists()) {
					dialog.setVisible(false);
					dialog.dispose();
					
					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame) SwingUtilities.getWindowAncestor(root), "Please wait...", "Creating Prosit CSV File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							LibraryFile library=new LibraryFile();
							library.openFile(libraryFile);
							PrositCSVWriter.writeCSV(library, defaultNCE, defaultCharge, false);
							return Nothing.NOTHING;
						}

						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify a library file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
		final JLabel text=new JLabel(AboutDialog.citationIcon);
		text.setText("<html><p style=\"font-size:10px; font-family: Helvetica, sans-serif\">"+"This function will extract out all peptides and charge states from an EncyclopeDIA library (.DLIB or .ELIB) and create an input file for Prosit. If you use this feature, please cite <a href=\"https://www.nature.com/articles/s41467-020-15346-1\">Searle et al, 2020</a>.");
		text.setBackground(Color.WHITE);
		text.setOpaque(true);
		text.addMouseListener(new MouseListener() {
			public void mouseClicked(MouseEvent e) {
			    try {
			         
			        Desktop.getDesktop().browse(new URI("https://www.nature.com/articles/s41467-020-15346-1"));
			         
			    } catch (IOException | URISyntaxException e1) {
			        e1.printStackTrace();
			    }
			}
			@Override
			public void mouseReleased(MouseEvent e) {
			}
			
			@Override
			public void mousePressed(MouseEvent e) {
			}
			
			@Override
			public void mouseExited(MouseEvent e) {
			}
			
			@Override
			public void mouseEntered(MouseEvent e) {
			}
		});

		mainpane.add(text, BorderLayout.NORTH);
		mainpane.add(options, BorderLayout.CENTER);
		mainpane.add(buttons, BorderLayout.SOUTH);
		mainpane.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10), BorderFactory.createTitledBorder("Parameters:")));
		
		dialog.getContentPane().add(mainpane, BorderLayout.CENTER);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		dialog.pack(); 
		dialog.setSize(500, 250);
		dialog.setVisible(true);
	}
	
	public static void convertLibraryForMS2PIP(Component root) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Convert Library to MS2PIP PEPREC", true);
		
		
		final FileChooserPanel libraryFileChooser=new FileChooserPanel(null, "Library", new SimpleFilenameFilter(".dlib", ".elib"), true);

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(libraryFileChooser);

		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				final File libraryFile=libraryFileChooser.getFile();
				if (libraryFile!=null&&libraryFile.exists()) {
					dialog.setVisible(false);
					dialog.dispose();
					
					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame) SwingUtilities.getWindowAncestor(root), "Please wait...", "Creating MS2PIP PEPREC File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							LibraryFile library=new LibraryFile();
							library.openFile(libraryFile);
							MS2PIPWriter.writeMS2PIP(library, false);
							return Nothing.NOTHING;
						}

						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify a library file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
		final JLabel text=new JLabel(AboutDialog.citationIcon);
		text.setText("<html><p style=\"font-size:10px; font-family: Helvetica, sans-serif\">"+"This function will extract out all peptides and charge states from an EncyclopeDIA library (.DLIB or .ELIB) and create an input file for MS2PIP. If you use this feature, please cite <a href=\"https://www.nature.com/articles/s41467-020-15346-1\">Searle et al, 2020</a>.");
		text.setBackground(Color.WHITE);
		text.setOpaque(true);
		text.addMouseListener(new MouseListener() {
			public void mouseClicked(MouseEvent e) {
			    try {
			         
			        Desktop.getDesktop().browse(new URI("https://www.nature.com/articles/s41467-020-15346-1"));
			         
			    } catch (IOException | URISyntaxException e1) {
			        e1.printStackTrace();
			    }
			}
			@Override
			public void mouseReleased(MouseEvent e) {
			}
			
			@Override
			public void mousePressed(MouseEvent e) {
			}
			
			@Override
			public void mouseExited(MouseEvent e) {
			}
			
			@Override
			public void mouseEntered(MouseEvent e) {
			}
		});

		mainpane.add(text, BorderLayout.NORTH);
		mainpane.add(options, BorderLayout.CENTER);
		mainpane.add(buttons, BorderLayout.SOUTH);
		mainpane.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10), BorderFactory.createTitledBorder("Parameters:")));
		
		dialog.getContentPane().add(mainpane, BorderLayout.CENTER);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		dialog.pack(); 
		dialog.setSize(500, 250);
		dialog.setVisible(true);
	}
	
	public static void convertSpectronaut(Component root, SearchParameters params) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Convert Prosit/Spectronaut CSV to Library", true);
		
		final FileChooserPanel csvFileChooser=new FileChooserPanel(null, "Spectronaut CSV/XLS", new SimpleFilenameFilter(".spectronaut", ".csv", ".tsv", ".txt", ".xls"), true);
		final FileChooserPanel fastaFileChooser=new FileChooserPanel(null, "FASTA", new SimpleFilenameFilter(".fas", ".fasta"), true);

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(csvFileChooser);
		options.add(fastaFileChooser);
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				final File csvFile=csvFileChooser.getFile();
				final File fastaFile=fastaFileChooser.getFile();
				
				if (csvFile!=null&&csvFile.exists()&&fastaFile!=null&&fastaFile.exists()) {
					dialog.setVisible(false);
					dialog.dispose();
					
					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame) SwingUtilities.getWindowAncestor(root), "Please wait...", "Reading Spectronaut CSV File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							SpectronautCSVToLibraryConverter.convertFromSpectronautCSV(csvFile, fastaFile, params);
							Logger.logLine("Finished reading "+csvFile.getName());
							return Nothing.NOTHING;
						}

						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify a Spectronaut CSV and a FASTA file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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

	public static void convertMaxQuantMSMSTXT(Component root, SearchParameters params) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Convert Maxquant msms.txt to Library", true);
		
		final FileChooserPanel csvFileChooser=new FileChooserPanel(null, "Maxquant msms.txt", new SimpleFilenameFilter("msms.txt"), true);
		final FileChooserPanel fastaFileChooser=new FileChooserPanel(null, "FASTA", new SimpleFilenameFilter(".fas", ".fasta"), true);

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(csvFileChooser);
		options.add(fastaFileChooser);
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				final File tsvFile=csvFileChooser.getFile();
				final File fastaFile=fastaFileChooser.getFile();
				
				if (tsvFile!=null&&tsvFile.exists()&&fastaFile!=null&&fastaFile.exists()) {
					dialog.setVisible(false);
					dialog.dispose();
					
					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame) SwingUtilities.getWindowAncestor(root), "Please wait...", "Reading Maxquant msms.txt File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							File libraryFile=new File(tsvFile.getAbsolutePath().substring(0, tsvFile.getAbsolutePath().lastIndexOf('.'))+LibraryFile.DLIB);
							MaxquantMSMSConverter.convertFromMSMSTSV(tsvFile, fastaFile, libraryFile, params);
							Logger.logLine("Finished reading "+tsvFile.getName());
							return Nothing.NOTHING;
						}

						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify a Maxquant msms.txt and a FASTA file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
	
	public static void convertMSP(Component root, SearchParameters params) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Convert NIST SPTXT/MSP to Library", true);
		
		final FileChooserPanel mspFileChooser=new FileChooserPanel(null, "SPTXT/MSP", new SimpleFilenameFilter(".msp", ".sptxt"), true);
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
				final File mspFile=mspFileChooser.getFile();
				final File fastaFile=fastaFileChooser.getFile();
				
				if (mspFile!=null&&mspFile.exists()&&fastaFile!=null&&fastaFile.exists()) {
					dialog.setVisible(false);
					dialog.dispose();
					
					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame) SwingUtilities.getWindowAncestor(root), "Please wait...", "Reading SPTXT/MSP File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							MSPReader.convertMSP(mspFile, fastaFile, params);
							Logger.logLine("Finished reading "+mspFile.getName());
							return Nothing.NOTHING;
						}

						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify a SPTXT/MSP and a FASTA file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
	
	public static void convertOpenSwathToELIB(Component root, SearchParameters params) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Convert OpenSwath TSV to Library", true);
		
		final FileChooserPanel tramlFileChooser=new FileChooserPanel(null, "OpenSwath TSV", new SimpleFilenameFilter(".tsv"), true);
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

				final File tramlFile=tramlFileChooser.getFile();
				final File fastaFile=fastaFileChooser.getFile();
				
				if (tramlFile!=null&&tramlFile.exists()&&fastaFile!=null&&fastaFile.exists()) {
					dialog.setVisible(false);
					dialog.dispose();
					
					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame) SwingUtilities.getWindowAncestor(root), "Please wait...", "Reading OpenSwath TSV File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							OpenSwathTSVToLibraryConverter.convertFromOpenSwathTSV(tramlFile, fastaFile, params);
							Logger.logLine("Finished reading "+tramlFile.getName());
							return Nothing.NOTHING;
						}

						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify an OpenSwath TSV and a FASTA file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
	
	public static void convertMS2PIPToELIB(Component root, SearchParameters params) {
		final JFrame frame = (JFrame)SwingUtilities.getRoot(root);
		final JDialog dialog=new JDialog(frame, "Convert MS2PIP to Library", true);

		final FileChooserPanel peprecFileChooser=new FileChooserPanel(null, "MS2PIP Input PEPREC", new SimpleFilenameFilter(".peprec"), true);
		final FileChooserPanel ms2pipCSVFileChooser=new FileChooserPanel(null, "MS2PIP Result CSV", new SimpleFilenameFilter(".csv"), true);
		final FileChooserPanel fastaFileChooser=new FileChooserPanel(null, "FASTA", new SimpleFilenameFilter(".fas", ".fasta"), true);

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(peprecFileChooser);
		options.add(ms2pipCSVFileChooser);
		options.add(fastaFileChooser);
		
		JPanel buttons=new JPanel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER));
		JButton okButton=new JButton("OK");
		okButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				final File peprecFile=peprecFileChooser.getFile();
				final File csvReportFile=ms2pipCSVFileChooser.getFile();
				final File fastaFile=fastaFileChooser.getFile();
				
				if (peprecFile!=null&&peprecFile.exists()&&csvReportFile!=null&&csvReportFile.exists()&&fastaFile!=null&&fastaFile.exists()) {
					dialog.setVisible(false);
					dialog.dispose();
					
					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame) SwingUtilities.getWindowAncestor(root), "Please wait...", "Reading MS2PIP File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							MS2PIPReader.convertMS2PIP(peprecFile, csvReportFile, fastaFile, params);
							Logger.logLine("Finished reading "+csvReportFile.getName());
							return Nothing.NOTHING;
						}

						@Override
						protected void doneForReal(Nothing t) {
						}
					};
					worker.execute();
				} else {
					JOptionPane.showMessageDialog(frame, "You must specify a MS2PIP PEPREC and CSV, and a FASTA file!", "Incomplete options!", JOptionPane.WARNING_MESSAGE, convertDBIcon);
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
		dialog.setSize(500, 210);
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

				final File tramlFile=tramlFileChooser.getFile();
				final File fastaFile=fastaFileChooser.getFile();
				
				if (tramlFile!=null&&tramlFile.exists()&&fastaFile!=null&&fastaFile.exists()) {
					if (tramlFile.length()>1024*1024*1000) {
					    int result = JOptionPane.showConfirmDialog(dialog, "This file is "+(tramlFile.length()/1024/1024)+" MB and will take a very long time to convert. Are you sure?",
					        "Warning: long conversion!", JOptionPane.OK_CANCEL_OPTION);
					    if (result==JOptionPane.CANCEL_OPTION) return;
					}
					
					dialog.setVisible(false);
					dialog.dispose();
					SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>((Frame) SwingUtilities.getWindowAncestor(root), "Please wait...", "Reading TraML File") {
						@Override
						protected Nothing doInBackgroundForReal() throws Exception {
							TraMLSAXToLibraryConverter.convertTraML(tramlFile, fastaFile, params);
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
		mainpane.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10), BorderFactory.createTitledBorder("Parameters: (ONLY USE THIS FOR SMALL LIBRARIES)")));
		
		dialog.getContentPane().add(mainpane, BorderLayout.CENTER);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		dialog.pack(); 
		dialog.setSize(500, 170);
		dialog.setVisible(true);
	}
}
