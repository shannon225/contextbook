package edu.washington.gs.maccoss.encyclopedia.gui.massspec;

import java.awt.Color;
import java.awt.Font;
import java.util.ArrayList;
import java.util.Optional;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.ui.RectangleInsets;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.ExtendedChartPanel;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;

public class ChromatogramCharter {
	public static ExtendedChartPanel createChart(Optional<ArrayList<XYTrace>> precursors, Optional<ArrayList<XYTrace>> fragments) {
		return createChart(precursors, fragments, 0.0, 0.0);
	}
	public static ExtendedChartPanel createChart(Optional<ArrayList<XYTrace>> precursors, Optional<ArrayList<XYTrace>> fragments, double globalMaxYPrecursor, double globalMaxYFragment) {
		ExtendedChartPanel fragmentPanel=Charter.getChart("Retention Time (min)", "MS2", false, globalMaxYFragment, 10, fragments.get().toArray(new XYTrace[0]));
		ExtendedChartPanel precursorPanel=Charter.getChart("Retention Time (min)", "MS1", false, globalMaxYPrecursor, 10, precursors.get().toArray(new XYTrace[0]));
		precursorPanel.getChart().getXYPlot().getRangeAxis().setInverted(true);

		ValueAxis domainAxis = fragmentPanel.getChart().getXYPlot().getDomainAxis();
		setFonts(domainAxis);
		setFonts(fragmentPanel.getChart().getXYPlot().getRangeAxis());
		setFonts(precursorPanel.getChart().getXYPlot().getRangeAxis());
		
		CombinedDomainXYPlot parent=new CombinedDomainXYPlot(domainAxis);
		parent.setGap(-1.0);
		parent.add(fragmentPanel.getChart().getXYPlot(), 3);
		parent.add(precursorPanel.getChart().getXYPlot(), 1);
		
		parent.setDomainGridlinesVisible(true);
		
		JFreeChart chart = new JFreeChart(parent);
		chart.setPadding(new RectangleInsets(10, 10, 10, 10));
		ExtendedChartPanel chartPanel=new ExtendedChartPanel(chart, false, fragmentPanel.getDivider());
		chartPanel.getChart().removeLegend();
		chartPanel.getChart().setBackgroundPaint(Color.white);
		chartPanel.setMinimumDrawWidth(0);
		chartPanel.setMinimumDrawHeight(0);
		chartPanel.setMaximumDrawWidth(Integer.MAX_VALUE);
		chartPanel.setMaximumDrawHeight(Integer.MAX_VALUE);
		
		return chartPanel;
	}

	private static void setFonts(ValueAxis axis) {
		axis.setTickLabelFont(new Font("News Gothic MT", Font.PLAIN, 12));
		axis.setLabelFont(new Font("News Gothic MT", Font.PLAIN, 12));
	}
	
	public static ArrayList<XYTrace> invert(ArrayList<XYTrace> traces) {
		ArrayList<XYTrace> newTraces=new ArrayList<>();
		for (XYTrace xyTrace : traces) {
			newTraces.add(xyTrace.rescaleY(-1f));
		}
		return newTraces;
	}
}
