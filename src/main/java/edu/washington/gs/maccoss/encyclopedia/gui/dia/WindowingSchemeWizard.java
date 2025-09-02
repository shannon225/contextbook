package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.BoxLayout;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.UIManager;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.datastructures.ScanRangeTracker;
import edu.washington.gs.maccoss.encyclopedia.datastructures.WindowSchemeGenerator;
import edu.washington.gs.maccoss.encyclopedia.gui.general.LabeledComponent;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;

public class WindowingSchemeWizard extends JPanel implements ActionListener, ChangeListener {
	private static final long serialVersionUID=1L;

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
				System.setProperty("com.apple.mrj.application.apple.menu.about.name", "Windowing Scheme Wizard");
				System.setProperty("apple.laf.useScreenMenuBar", "true");
				break;

			default:
				break;
		}

		EventQueue.invokeLater(new Runnable() {
			@Override
			public void run() {
				final JFrame f=new JFrame("Windowing Scheme Wizard");
				f.addWindowListener(new WindowAdapter() {
					public void windowClosing(WindowEvent e) {
						System.exit(0);
					}
				});

				f.getContentPane().add(new WindowingSchemeWizard(), BorderLayout.CENTER);

				f.pack();
				f.setSize(new Dimension(1900, 1030)); // for 1920x1080
				f.setVisible(true);
			}
		});

		Logger.logLine("Launching Windowing Scheme Wizard");
	}

	// follows indexing in WindowSchemeGenerator
	private static final String[] WINDOWING_SCHEME_ITEMS=new String[] {"Normal DIA", "Staggered/Overlap DIA", "Variable Width DIA"};
	private static JCheckBox isPhospho=new JCheckBox();
	private final SpinnerModel numberOfWindows=new SpinnerNumberModel(25, 1, 500, 1);
	private final SpinnerModel startMz=new SpinnerNumberModel(400, 1, 2000, 1);
	private final SpinnerModel stopMz=new SpinnerNumberModel(1000, 1, 2000, 1);
	private final SpinnerModel marginWidth=new SpinnerNumberModel(0, 0, 5, 0.1);
	private final JComboBox<String> windowingScheme=new JComboBox<String>(WINDOWING_SCHEME_ITEMS);
	private final JSplitPane structureSplit=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	private final JSplitPane tableSplit=new JSplitPane(JSplitPane.VERTICAL_SPLIT);
	private final JTable table;
	private final ScanRangeTableModel model;
	private final TableRowSorter<TableModel> rowSorter;
	
	public WindowingSchemeWizard() {
		super(new BorderLayout());

		JPanel options=new JPanel();
		options.setLayout(new BoxLayout(options, BoxLayout.PAGE_AXIS));
		options.add(new LabeledComponent("<p style=\"font-size:12px; font-family: Helvetica, sans-serif\"><b>Parameters", new JLabel()));

		JSpinner c1=new JSpinner(numberOfWindows);
		JSpinner c2=new JSpinner(startMz);
		JSpinner c3=new JSpinner(stopMz);
		JSpinner c4=new JSpinner(marginWidth);
		
		options.add(new LabeledComponent("Windowing Scheme", windowingScheme));
		options.add(new LabeledComponent("Phospho Enriched", isPhospho));
		options.add(new LabeledComponent("Number Of Windows", c1));
		options.add(new LabeledComponent("Start m/z", c2));
		options.add(new LabeledComponent("Stop m/z", c3));
		options.add(new LabeledComponent("Margin Width", c4));

		windowingScheme.addActionListener(this);
		isPhospho.addActionListener(this);
		c1.addChangeListener(this);
		c2.addChangeListener(this);
		c3.addChangeListener(this);
		c4.addChangeListener(this);

		model=new ScanRangeTableModel();
		table=new JTable(model);
		
		rowSorter=new TableRowSorter<TableModel>(table.getModel());
		table.setRowSorter(rowSorter);
		
		tableSplit.setTopComponent(options);
		tableSplit.setBottomComponent(new JScrollPane(table));
		
		structureSplit.setLeftComponent(tableSplit);
		structureSplit.setRightComponent(new JPanel());
		this.add(structureSplit, BorderLayout.CENTER);
		
		updateStructure();
	}

	@Override
	public void stateChanged(ChangeEvent e) {
		updateStructure();
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		updateStructure();
	}
	
	private void updateStructure() {
		int structureLocation=structureSplit.getDividerLocation();
		if (structureLocation<=5) {
			structureLocation=400;
		}
		int tableLocation=tableSplit.getDividerLocation();
		if (tableLocation<=5) {
			tableLocation=200;
		}
		
		ScanRangeTracker scanRange=generateScanRangeTracker(isPhospho.isSelected());
		ChartPanel structureChart=MzmlStructureCharter.getStructureChart(scanRange, true);
		structureChart.restoreAutoBounds();
		structureSplit.setRightComponent(structureChart);
		model.updateEntries(scanRange);

		structureSplit.setDividerLocation(structureLocation);
		tableSplit.setDividerLocation(tableLocation);
	}
	
	public ScanRangeTracker generateScanRangeTracker(boolean isPhospho) {
		int start=((Number)startMz.getValue()).intValue();
		int stop=((Number)stopMz.getValue()).intValue();
		int numWindows=((Number)numberOfWindows.getValue()).intValue();
		float margin=((Number)marginWidth.getValue()).floatValue();
		int windowingSchemeIndex=windowingScheme.getSelectedIndex();
		if (isPhospho&&windowingSchemeIndex==WindowSchemeGenerator.VARIABLE_WIDTH_DIA) {
			windowingSchemeIndex=WindowSchemeGenerator.VARIABLE_WIDTH_PHOSPHO_DIA;
		}

		return WindowSchemeGenerator.generateWindowingScheme(windowingSchemeIndex, start, stop, numWindows, margin, isPhospho);
	}
}
