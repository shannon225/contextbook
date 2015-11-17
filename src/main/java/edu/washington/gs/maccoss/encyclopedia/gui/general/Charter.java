package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.geom.Ellipse2D;
import java.security.InvalidParameterException;

import javax.swing.JFrame;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.AbstractXYItemRenderer;
import org.jfree.chart.renderer.xy.XYAreaRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;

public class Charter {
	public static void main(String[] args) {
		XYTrace trace=new XYTrace(new double[] {114.0913405, 147.1128042, 227.1754045, 244.1655682, 355.2339825, 359.1925112, 458.2609252, 484.2765755, 515.2823892, 541.2980395, 640.3664535,
				644.3249822, 755.3933965, 772.3835602, 852.4461605, 885.4676242, 980.5411235, 998.5516882}, new double[] {0.021099463, 0.00721319, 0.10845732, 0.116413645, 0.39157316, 0.1849763,
				0.443399, 0.35894206, 0.43697295, 0.47858942, 0.5025189, 0.34656474, 0.26218376, 0.27163184, 0.2108471, 0.23929471, 0.12108889, 0.12206937}, GraphType.spectrum, "Trace");
		XYTrace trace2=new XYTrace(new double[] {114.0913405, 147.1128042, 227.1754045, 244.1655682, 355.2339825, 359.1925112, 458.2609252, 484.2765755, 515.2823892, 541.2980395, 640.3664535,
				644.3249822, 755.3933965, 772.3835602, 852.4461605, 885.4676242, 980.5411235, 998.5516882}, new double[] {0.021099463, 0.00721319, 0.10845732, 0.116413645, 0.39157316, 0.1849763,
				0.443399, 0.35894206, 0.43697295, 0.47858942, 0.5025189, 0.34656474, 0.26218376, 0.27163184, 0.2108471, 0.23929471, 0.12108889, 0.12206937}, GraphType.line, "Trace2");
		
		launchChart("M/Z", "Intensity", trace, trace2);
	}

	public static void launchChart(String xAxis, String yAxis, XYTrace... traces) {
		final JFrame f=new JFrame("Runs vs Reference");
		f.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});

		f.getContentPane().add(getChart(xAxis, yAxis, traces), BorderLayout.CENTER);

		f.pack();
		f.setSize(new Dimension(1000, 770));
		f.setVisible(true);
	}

	public static ChartPanel getChart(String xAxis, String yAxis, XYTrace... traces) {
		NumberAxis numberaxis=new NumberAxis(xAxis);
		numberaxis.setAutoRangeIncludesZero(false);
		NumberAxis numberaxis1=new NumberAxis(yAxis);
		numberaxis1.setAutoRangeIncludesZero(false);

		XYPlot plot=new XYPlot();
		plot.setDomainAxis(numberaxis);
		plot.setRangeAxis(numberaxis1);

		int count=0;
		for (XYTrace trace : traces) {
			AbstractXYItemRenderer renderer=new XYLineAndShapeRenderer();
			switch (trace.getType()) {
				case area:
					renderer=new XYAreaRenderer();
					renderer.setSeriesStroke(0, new BasicStroke(2.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
					break;

				case line:
					renderer=new XYLineAndShapeRenderer();
					renderer.setSeriesStroke(0, new BasicStroke(2.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
					((XYLineAndShapeRenderer)renderer).setBaseShapesVisible(false);

					break;

				case point:
					renderer=new XYLineAndShapeRenderer();
					renderer.setSeriesShape(0, new Ellipse2D.Double(0, 0, 5, 5));
					((XYLineAndShapeRenderer)renderer).setBaseLinesVisible(false);

					break;

				case spectrum:
					renderer=new XYLineAndShapeRenderer();
					((XYLineAndShapeRenderer)renderer).setBaseShapesVisible(false);
					renderer.setBasePaint(Color.black);

					break;

				default:
					throw new InvalidParameterException("unsupported graphing type!");
			}

			Pair<double[], double[]> values=trace.toArrays();
			double[] x=values.x;
			double[] y=values.y;
			XYSeriesCollection dataset=new XYSeriesCollection();
			switch (trace.getType()) {
				case area:
				case line:
				case point:
					XYSeries series=new XYSeries(trace.getName());
					for (int i=0; i<x.length; i++) {
						series.add(x[i], y[i]);
					}
					dataset.addSeries(series);
					break;

				case spectrum:
					for (int i=0; i<x.length; i++) {
						XYSeries peakSeries=new XYSeries(x[i]);
						peakSeries.add(x[i], 0);
						peakSeries.add(x[i], y[i]);
						dataset.addSeries(peakSeries);
						renderer.setSeriesStroke(i, new BasicStroke(2.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
						renderer.setSeriesPaint(i, new Color(26, 148, 49));
					}
					break;

				default:
					throw new InvalidParameterException("unsupported graphing type!");
			}

			plot.setDataset(count, dataset);
			plot.setRenderer(count, renderer);

			count++;
		}

		plot.setBackgroundPaint(Color.white);
		plot.setDomainGridlinePaint(Color.gray);
		plot.setRangeGridlinePaint(Color.gray);
		JFreeChart chart=new JFreeChart(plot);
		chart.setBackgroundPaint(Color.white);
		ChartPanel chartPanel=new ChartPanel(chart, false);
		chartPanel.getChart().removeLegend();

		NumberAxis rangeAxis=(NumberAxis)((XYPlot)plot).getRangeAxis();

		Font font=new Font("News Gothic MT", Font.PLAIN, 24);
		Font font2=new Font("News Gothic MT", Font.PLAIN, 32);
		rangeAxis.setLabelFont(font2);
		rangeAxis.setTickLabelFont(font);

		NumberAxis domainAxis=(NumberAxis)((XYPlot)plot).getDomainAxis();
		if (domainAxis!=null) {
			domainAxis.setLabelFont(font2);
			domainAxis.setTickLabelFont(font);
		}

		return chartPanel;
	}

}
