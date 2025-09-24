package edu.washington.gs.maccoss.encyclopedia.gui.massspec;

import java.awt.Color;
import java.awt.Font;
import java.util.ArrayList;
import java.util.Optional;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.ui.RectangleInsets;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.ExtendedChartPanel;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;

public class ChromatogramCharter {
	private static final int DEFAULT_FONT_SIZE = 16;

	public static ExtendedChartPanel createChart(Optional<ArrayList<XYTrace>> precursors, Optional<ArrayList<XYTrace>> fragments) {
		return createChart(precursors, fragments, 0.0, 0.0, DEFAULT_FONT_SIZE);
	}
	public static ExtendedChartPanel createChart(Optional<ArrayList<XYTrace>> precursors, Optional<ArrayList<XYTrace>> fragments, int fontSize) {
		return createChart(precursors, fragments, 0.0, 0.0, fontSize);
	}

	public static ExtendedChartPanel createChart(Optional<ArrayList<XYTrace>> precursors, Optional<ArrayList<XYTrace>> fragments, double globalMaxYPrecursor, double globalMaxYFragment) {
		return createChart(precursors, fragments, globalMaxYPrecursor, globalMaxYFragment, DEFAULT_FONT_SIZE);
	}
	public static ExtendedChartPanel createChart(Optional<ArrayList<XYTrace>> precursors, Optional<ArrayList<XYTrace>> fragments, double globalMaxYPrecursor, double globalMaxYFragment, int fontSize) {
		if (!precursors.isPresent()&&!fragments.isPresent()) {
			throw new EncyclopediaException("Precursors and fragments can't both be missing!");
		} else if (!precursors.isPresent()) {
			ExtendedChartPanel fragmentPanel=Charter.getChart("Retention Time (min)", "MS2", false, globalMaxYFragment, fontSize, fragments.get().toArray(new XYTrace[0]));
			return fragmentPanel;
		} else if (!fragments.isPresent()) {
			ExtendedChartPanel precursorPanel=Charter.getChart("Retention Time (min)", "MS1", false, globalMaxYPrecursor, fontSize, precursors.get().toArray(new XYTrace[0]));
			return precursorPanel;
		}
		ExtendedChartPanel fragmentPanel=Charter.getChart("Retention Time (min)", "MS2", false, globalMaxYFragment, fontSize, fragments.get().toArray(new XYTrace[0]));
		ExtendedChartPanel precursorPanel=Charter.getChart("Retention Time (min)", "MS1", false, globalMaxYPrecursor, fontSize, precursors.get().toArray(new XYTrace[0]));
		precursorPanel.getChart().getXYPlot().getRangeAxis().setInverted(true);

		ValueAxis domainAxis = fragmentPanel.getChart().getXYPlot().getDomainAxis();
		setFonts(domainAxis, fontSize);
		setFonts(fragmentPanel.getChart().getXYPlot().getRangeAxis(), fontSize);
		setFonts(precursorPanel.getChart().getXYPlot().getRangeAxis(), fontSize);
		
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
		
		ArrayList<XYTrace> allTraces=new ArrayList<XYTrace>();
		if (precursors.isPresent()) allTraces.addAll(precursors.get());
		if (fragments.isPresent()) allTraces.addAll(fragments.get());
		
		Charter.addCopyDataMenu("Retention Time (min)", chartPanel, allTraces.toArray(new XYTrace[0]));
		
		return chartPanel;
	}

	private static void setFonts(ValueAxis axis, int fontSize) {
		axis.setTickLabelFont(new Font(Charter.BASE_FONT_NAME, Font.PLAIN, fontSize));
		axis.setLabelFont(new Font(Charter.BASE_FONT_NAME, Font.PLAIN, fontSize));
	}
	
	public static ArrayList<XYTrace> invert(ArrayList<XYTrace> traces) {
		ArrayList<XYTrace> newTraces=new ArrayList<>();
		for (XYTrace xyTrace : traces) {
			newTraces.add(xyTrace.rescaleY(-1f));
		}
		return newTraces;
	}
}
