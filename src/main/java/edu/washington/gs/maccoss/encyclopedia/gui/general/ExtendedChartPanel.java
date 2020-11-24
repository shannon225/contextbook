package edu.washington.gs.maccoss.encyclopedia.gui.general;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;

public class ExtendedChartPanel extends ChartPanel {
	private static final long serialVersionUID = 1L;
	
	private final double divider;

	public ExtendedChartPanel(JFreeChart chart, boolean useBuffer, double divider) {
		super(chart, useBuffer);
		this.divider=divider;
	}

	public ExtendedChartPanel(JFreeChart chart, double divider) {
		super(chart);
		this.divider=divider;
	}

	public double getDivider() {
		return divider;
	}
}
