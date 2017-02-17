package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.FileOutputStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Map;
import java.util.Map.Entry;

import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JTabbedPane;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.CategoryLabelPositions;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.PaintScale;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.chart.renderer.category.StandardBarPainter;
import org.jfree.chart.renderer.xy.AbstractXYItemRenderer;
import org.jfree.chart.renderer.xy.XYAreaRenderer;
import org.jfree.chart.renderer.xy.XYBlockRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.TextAnchor;

import com.itextpdf.awt.PdfGraphics2D;
import com.itextpdf.text.Document;
import com.itextpdf.text.Rectangle;
import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfTemplate;
import com.itextpdf.text.pdf.PdfWriter;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYZTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution;

public class Charter {
	public static void main(String[] args) {
		XYTraceInterface trace=new XYTrace(
				new double[] { 114.0913405, 147.1128042, 227.1754045, 244.1655682, 355.2339825, 359.1925112, 458.2609252, 484.2765755, 515.2823892, 541.2980395, 640.3664535, 644.3249822, 755.3933965,
						772.3835602, 852.4461605, 885.4676242, 980.5411235, 998.5516882 },
				new double[] { 0.021099463, 0.00721319, 0.10845732, 0.116413645, 0.39157316, 0.1849763, 0.443399, 0.35894206, 0.43697295, 0.47858942, 0.5025189, 0.34656474, 0.26218376, 0.27163184,
						0.2108471, 0.23929471, 0.12108889, 0.12206937 },
				GraphType.spectrum, "Trace");
		XYTraceInterface trace2=new XYTrace(
				new double[] { 114.0913405, 147.1128042, 227.1754045, 244.1655682, 355.2339825, 359.1925112, 458.2609252, 484.2765755, 515.2823892, 541.2980395, 640.3664535, 644.3249822, 755.3933965,
						772.3835602, 852.4461605, 885.4676242, 980.5411235, 998.5516882 },
				new double[] { 0.021099463, 0.00721319, 0.10845732, 0.116413645, 0.39157316, 0.1849763, 0.443399, 0.35894206, 0.43697295, 0.47858942, 0.5025189, 0.34656474, 0.26218376, 0.27163184,
						0.2108471, 0.23929471, 0.12108889, 0.12206937 },
				GraphType.line, "Trace2");
		XYTraceInterface trace3=new XYTrace(new double[] { 114.0913405, 147.1128042, 227.1754045, 244.1655682, 355.2339825, 359.1925112, 458.2609252, 484.2765755, 515.2823892, 541.2980395, 640.3664535,
				644.3249822, 755.3933965, 772.3835602, 852.4461605, 885.4676242, 980.5411235, 998.5516882 }, new double[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, GraphType.spectrum,
				"Trace3");

		XYTraceInterface trace4=new XYTrace(
				new double[] {347.85, 426.35, 551.4, 552.13, 553.31, 553.87, 569.2, 571.3, 621.01, 621.88, 622.21, 637.77, 638.31, 639.15, 640.14, 648.13, 685.4, 703.51, 703.92, 704.61, 705.83,
						706.33, 735.48, 735.9, 736.43, 755.31, 756.12, 757.03, 770.11, 772.55, 772.93, 774.72, 800.35, 824.93, 825.41, 825.85, 826.13, 826.86, 827.46, 843.2, 852.26, 852.86, 853.71,
						854.45, 870.63, 872.11, 912.62},
				new double[] {1713, 1166, 1941, 1890, 1873, 1465, 1151, 1456, 1936, 2725, 1253, 1352, 2858, 1810, 1334, 2466, 1279, 9072, 4986, 5070, 1320, 1220, 1202, 1550, 5121, 2735, 1452, 1084,
						1136, 1661, 1544, 1620, 1210, 1770, 1297, 2086, 1635, 1117, 1432, 1652, 2143, 3766, 5043, 2124, 54156, 1310, 1240},
				GraphType.dashedline, "Trace4");

		ChartPanel chart=getChart("M/Z", "Intensity", false, trace, trace2, trace3);
		chart=getChart("M/Z", "Intensity", false, trace4);
		launchChart(chart, "Title!");

		// writeAsPDF(chart.getChart(), new
		// File("/Users/searleb/Documents/projects/encyclopedia/mzml/test.pdf"),
		// new Dimension(792, 612));
	}

	public static void launchComponent(JComponent comp, String title, Dimension dim) {
		final JFrame f=new JFrame(title);
		f.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});

		f.getContentPane().add(comp, BorderLayout.CENTER);

		f.pack();
		f.setSize(dim);
		f.setVisible(true);
	}

	public static void launchChart(ChartPanel chart, String title) {
		launchComponent(chart, title, new Dimension(792, 612));
	}

	public static void launchChart(Spectrum trace) {
		launchChart(trace, trace.getSpectrumName());
	}

	public static void launchChart(Spectrum trace, String title) {
		launchComponent(getChart(trace), title, new Dimension(1000, 500));
	}

	public static void launchCharts(String title, Map<String, ChartPanel> panelMap) {
		launchComponent(getTabbedChartPane(panelMap), title, new Dimension(792, 612));
	}

	public static JTabbedPane getTabbedChartPane(Map<String, ChartPanel> panelMap) {
		JTabbedPane tabs=new JTabbedPane();
		for (Entry<String, ChartPanel> entry : panelMap.entrySet()) {
			tabs.addTab(entry.getKey(), entry.getValue());
		}
		return tabs;
	}

	public static void launchChart(String xAxis, String yAxis, boolean displayLegend, Dimension dim, XYTraceInterface... traces) {
		launchComponent(getChart(xAxis, yAxis, displayLegend, traces), xAxis, dim);
	}

	public static void launchChart(String xAxis, String yAxis, boolean displayLegend, XYTraceInterface... traces) {
		launchComponent(getChart(xAxis, yAxis, displayLegend, traces), xAxis, new Dimension(792, 612));
	}

	public static void launchChart(String xAxis, String yAxis, boolean displayLegend, XYZTrace dataset) {
		launchComponent(getChart(xAxis, yAxis, displayLegend, dataset), xAxis, new Dimension(792, 612));
	}

	public static void writeAsPDF(File f, String xAxis, String yAxis, boolean displayLegend, XYTraceInterface... traces) {
		//Dimension d=new Dimension(792, 612);
		Dimension d=new Dimension(400, 300);
		
		writeAsPDF(getChart(xAxis, yAxis, displayLegend, traces).getChart(), f, d);
	}

	public static void writeAsPDF(JFreeChart chart, File f, Dimension d) {
		try {
			Rectangle pagesize=new Rectangle(d.width, d.height);
			Document document=new Document(pagesize);
			FileOutputStream os=new FileOutputStream(f);
			PdfWriter writer=PdfWriter.getInstance(document, os);
			document.open();
			PdfContentByte canvas=writer.getDirectContent();
			PdfTemplate template=canvas.createTemplate(d.width, d.height);
			Graphics2D g2d=new PdfGraphics2D(template, d.width, d.height);

			Rectangle2D r2D=new Rectangle2D.Double(0, 0, d.width, d.height);
			chart.draw(g2d, r2D);
			g2d.dispose();
			canvas.addTemplate(template, 0, 0);
			document.close();
			os.close();
		} catch (Exception e) {
			Logger.errorException(e);
		}
	}

	public static ChartPanel getChart(String xAxisName, String yAxisName, boolean displayLegend, XYZTrace dataset) {
		NumberAxis xAxis=new NumberAxis(xAxisName);
		xAxis.setAutoRangeIncludesZero(false);
		NumberAxis yAxis=new NumberAxis(yAxisName);
		yAxis.setAutoRangeIncludesZero(false);

		XYBlockRenderer renderer=new XYBlockRenderer();
		renderer.setBlockHeight(5);
		renderer.setBlockWidth(5);
		PaintScale scale=new PeakPaintScale(dataset.getMinZ(), dataset.getMaxZ());
		renderer.setPaintScale(scale);

		XYPlot plot=new XYPlot(dataset, xAxis, yAxis, renderer);
		plot.setBackgroundPaint(Color.white);
		plot.setDomainGridlinesVisible(false);
		plot.setRangeGridlinePaint(Color.white);
		plot.setRangeGridlinesVisible(false);
		JFreeChart chart=new JFreeChart(plot);
		chart.removeLegend();
		chart.setBackgroundPaint(Color.white);

		NumberAxis rangeAxis=(NumberAxis) plot.getRangeAxis();

		Font font=new Font("News Gothic MT", Font.PLAIN, 24);
		Font font2=new Font("News Gothic MT", Font.PLAIN, 32);
		Font font3=new Font("News Gothic MT", Font.PLAIN, 18);
		rangeAxis.setLabelFont(font2);
		rangeAxis.setTickLabelFont(font);

		ChartPanel chartPanel=new ChartPanel(chart, false);
		if (!displayLegend) {
			chartPanel.getChart().removeLegend();
		} else {
			chartPanel.getChart().getLegend().setItemFont(font3);
		}

		NumberAxis domainAxis=(NumberAxis) ((XYPlot) plot).getDomainAxis();
		if (domainAxis!=null) {
			domainAxis.setLabelFont(font2);
			domainAxis.setTickLabelFont(font);
		}
		return chartPanel;
	}

	public static ChartPanel getChart(LibraryEntry trace) {
		ChartPanel chart=getChart("M/Z", "Intensity", false, trace);
		chart.getChart().setTitle(trace.getSpectrumName());
		return chart;
	}

	public static ChartPanel getChart(Spectrum trace) {
		ChartPanel chart=getChart("M/Z", "Intensity", false, new XYTrace(trace));
		chart.getChart().setTitle(trace.getSpectrumName());
		return chart;
	}

	public static ChartPanel getChart(Distribution dist, Range range) {
		int n=100;

		ArrayList<XYPoint> points=new ArrayList<XYPoint>();
		for (int i=0; i<n; i++) {
			float value=(i/(float) n)*range.getRange()+range.getStart();
			points.add(new XYPoint(value, dist.getProbability(value)));
		}

		return getChart("Value", "Probability", false, new XYTrace(points, GraphType.line, dist.getName()));
	}
	private static final DecimalFormat MASS_FORMAT = new DecimalFormat(".#");
	
	{
		BarRenderer.setDefaultBarPainter(new StandardBarPainter());
		BarRenderer.setDefaultShadowsVisible(false);
	}
	public static ChartPanel getBarChart(String title, String xAxis, String yAxis, String[] categories, float[] values) {
		assert (categories.length==values.length);
		boolean displayLegend=false;

		DefaultCategoryDataset dataset=new DefaultCategoryDataset();
		for (int i=0; i<values.length; i++) {
			dataset.addValue(values[i], xAxis, categories[i]);
		}
		JFreeChart barChart=ChartFactory.createBarChart(title, xAxis, yAxis, dataset, PlotOrientation.VERTICAL, false, true, false);

		CategoryPlot plot=barChart.getCategoryPlot();
	    ((BarRenderer)plot.getRenderer()).setBarPainter(new StandardBarPainter());

		plot.setBackgroundPaint(Color.white);
		plot.setDomainGridlinePaint(Color.white);//gray);
		plot.setDomainGridlinesVisible(false);
		plot.setRangeGridlinePaint(Color.white);//gray);
		plot.setRangeGridlinesVisible(false);
		JFreeChart chart=new JFreeChart(plot);
		chart.setBackgroundPaint(Color.white);

		NumberAxis rangeAxis=(NumberAxis) ((CategoryPlot) plot).getRangeAxis();

		Font font=new Font("News Gothic MT", Font.PLAIN, 24);
		Font font2=new Font("News Gothic MT", Font.PLAIN, 32);
		Font font3=new Font("News Gothic MT", Font.PLAIN, 18);
		font=new Font("News Gothic MT", Font.PLAIN, 12);
		font2=new Font("News Gothic MT", Font.PLAIN, 16);
		font3=new Font("News Gothic MT", Font.PLAIN, 16);
		rangeAxis.setLabelFont(font2);
		rangeAxis.setTickLabelFont(font);

		CategoryAxis domainAxis=(CategoryAxis)((CategoryPlot)plot).getDomainAxis();
		domainAxis.setCategoryLabelPositions(CategoryLabelPositions.UP_90);
	    domainAxis.setMaximumCategoryLabelLines(1);
	    domainAxis.setMaximumCategoryLabelWidthRatio(2f);
		if (domainAxis!=null) {
			domainAxis.setLabelFont(font);
			domainAxis.setTickLabelFont(font);
		}

		ChartPanel chartPanel=new ChartPanel(chart, false);
		if (!displayLegend) {
			chartPanel.getChart().removeLegend();
		} else {
			chartPanel.getChart().getLegend().setItemFont(font3);
		}

		chartPanel.setMinimumDrawWidth(0);
		chartPanel.setMinimumDrawHeight(0);
		chartPanel.setMaximumDrawWidth(Integer.MAX_VALUE);
		chartPanel.setMaximumDrawHeight(Integer.MAX_VALUE);
		return chartPanel;
	}

	public static ChartPanel getChart(String xAxis, String yAxis, boolean displayLegend, XYTraceInterface... traces) {
		NumberAxis numberaxis=new NumberAxis(xAxis);
		numberaxis.setAutoRangeIncludesZero(false);
		NumberAxis numberaxis1=new NumberAxis(yAxis);
		numberaxis1.setAutoRangeIncludesZero(false);

		XYPlot plot=new XYPlot();
		plot.setDomainAxis(numberaxis);
		plot.setRangeAxis(numberaxis1);

		int count=0;
		for (XYTraceInterface trace : traces) {
			AbstractXYItemRenderer renderer=new XYLineAndShapeRenderer();
			switch (trace.getType()) {
			case area:
				renderer=new XYAreaRenderer();
				renderer.setSeriesStroke(0, new BasicStroke(2.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
				break;

			case line:
				renderer=new XYLineAndShapeRenderer();
				renderer.setSeriesStroke(0, new BasicStroke(trace.getThickness().orElse(2.0f), BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
				((XYLineAndShapeRenderer) renderer).setBaseShapesVisible(false);

				break;

			case dashedline:
				renderer=new XYLineAndShapeRenderer();
				renderer.setSeriesStroke(0, new BasicStroke(trace.getThickness().orElse(2.0f), BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 0.0f, new float[] {3.0f, 5.0f}, 0.0f));
				((XYLineAndShapeRenderer) renderer).setDrawSeriesLineAsPath(true);
				((XYLineAndShapeRenderer) renderer).setBaseShapesVisible(false);

				break;

			case point:
				renderer=new XYLineAndShapeRenderer();
				renderer.setSeriesShape(0, new Ellipse2D.Double(0, 0, 3, 3));
				((XYLineAndShapeRenderer) renderer).setBaseLinesVisible(false);

				break;

			case tinypoint:
				renderer=new XYLineAndShapeRenderer();
				renderer.setSeriesShape(0, new Ellipse2D.Double(0, 0, 1, 1));
				((XYLineAndShapeRenderer) renderer).setBaseLinesVisible(false);

				break;

			case spectrum:
				renderer=new XYLineAndShapeRenderer();
				((XYLineAndShapeRenderer) renderer).setBaseShapesVisible(false);
				renderer.setBasePaint(Color.black);

				break;

			default:
				throw new EncyclopediaException("unsupported graphing type!");
			}

			Pair<double[], double[]> values=trace.toArrays();
			double[] x=values.x;
			double[] y=values.y;
			XYSeriesCollection dataset=new XYSeriesCollection();
			switch (trace.getType()) {
			case area:
			case line:
			case dashedline:
			case point:
			case tinypoint:
				XYSeries series=new XYSeries(trace.getName());
				for (int i=0; i<x.length; i++) {
					if (!Double.isNaN(x[i])&&!Double.isNaN(y[i])) {
						series.add(x[i], y[i]);
					}
				}
				dataset.addSeries(series);
				break;

			case spectrum:
				double yThreshold=General.max(y)*0.2;
				if (trace instanceof AnnotatedLibraryEntry) {
					AnnotatedLibraryEntry entry=(AnnotatedLibraryEntry)(Spectrum)trace;
					FragmentIon[] annotations=entry.getIonAnnotations();
					
					for (int i=0; i<x.length; i++) {
						if (!Double.isNaN(x[i])&&!Double.isNaN(y[i])) {
							XYSeries peakSeries=new XYSeries(x[i]);
							peakSeries.add(x[i], 0);
							peakSeries.add(x[i], y[i]);
							dataset.addSeries(peakSeries);

							if (annotations[i]!=null) {
								renderer.setSeriesStroke(i, IonType.getStroke(annotations[i].type));
								renderer.setSeriesPaint(i, IonType.getColor(annotations[i].type));
								
								XYTextAnnotation xytextannotation = new XYTextAnnotation(annotations[i].toString(), x[i], y[i]);
								xytextannotation.setPaint(IonType.getColor(annotations[i].type));
						        xytextannotation.setFont(IonType.getFont(annotations[i].type));
						        xytextannotation.setTextAnchor(TextAnchor.BOTTOM_CENTER);
						        plot.addAnnotation(xytextannotation);
							} else {
								renderer.setSeriesStroke(i, IonType.missingStroke);
								renderer.setSeriesPaint(i, IonType.missingColor);
								
								if (y[i]>yThreshold) {
									XYTextAnnotation xytextannotation = new XYTextAnnotation(MASS_FORMAT.format(x[i]), x[i], y[i]);
									xytextannotation.setPaint(IonType.missingColor);
							        xytextannotation.setFont(IonType.missingAnnotationFont);
							        xytextannotation.setTextAnchor(TextAnchor.BOTTOM_CENTER);
							        plot.addAnnotation(xytextannotation);
								}
							}
						}
					}
				} else {
					for (int i=0; i<x.length; i++) {
						if (!Double.isNaN(x[i])&&!Double.isNaN(y[i])) {
							XYSeries peakSeries=new XYSeries(x[i]);
							peakSeries.add(x[i], 0);
							peakSeries.add(x[i], y[i]);
							dataset.addSeries(peakSeries);
							renderer.setSeriesStroke(i, IonType.missingStroke);
							renderer.setSeriesPaint(i, IonType.missingColor);
							
							if (y[i]>yThreshold) {
								XYTextAnnotation xytextannotation = new XYTextAnnotation(MASS_FORMAT.format(x[i]), x[i], y[i]);
								xytextannotation.setPaint(IonType.missingColor);
						        xytextannotation.setFont(IonType.missingAnnotationFont);
						        xytextannotation.setTextAnchor(TextAnchor.BOTTOM_CENTER);
						        plot.addAnnotation(xytextannotation);
							}
						}
					}
				}
				break;

			default:
				throw new EncyclopediaException("unsupported graphing type!");
			}
			if (trace.getColor().isPresent()) {
				renderer.setSeriesPaint(0, trace.getColor().get());
				renderer.setBasePaint(trace.getColor().get());
			}

			plot.setDataset(count, dataset);
			plot.setRenderer(count, renderer);

			count++;
		}

		plot.setBackgroundPaint(Color.white);
		plot.setDomainGridlinePaint(Color.white);//gray);
		plot.setDomainGridlinesVisible(false);
		plot.setRangeGridlinePaint(Color.white);//gray);
		plot.setRangeGridlinesVisible(false);
		JFreeChart chart=new JFreeChart(plot);
		chart.setBackgroundPaint(Color.white);

		NumberAxis rangeAxis=(NumberAxis) ((XYPlot) plot).getRangeAxis();

		Font font=new Font("News Gothic MT", Font.PLAIN, 24);
		Font font2=new Font("News Gothic MT", Font.PLAIN, 32);
		Font font3=new Font("News Gothic MT", Font.PLAIN, 18);
		font=new Font("News Gothic MT", Font.PLAIN, 16);
		font2=new Font("News Gothic MT", Font.PLAIN, 16);
		font3=new Font("News Gothic MT", Font.PLAIN, 16);
		rangeAxis.setLabelFont(font2);
		rangeAxis.setTickLabelFont(font);

		NumberAxis domainAxis=(NumberAxis) ((XYPlot) plot).getDomainAxis();
		if (domainAxis!=null) {
			domainAxis.setLabelFont(font2);
			domainAxis.setTickLabelFont(font);
		}

		ChartPanel chartPanel=new ChartPanel(chart, false);
		if (!displayLegend) {
			chartPanel.getChart().removeLegend();
		} else {
			chartPanel.getChart().getLegend().setItemFont(font3);
		}
		//rangeAxis.setTickUnit(new NumberTickUnit(20)); 
		//domainAxis.setTickUnit(new NumberTickUnit(20));
		//rangeAxis.setRange(0, 135);
		//domainAxis.setRange(15, 135);

		chartPanel.setMinimumDrawWidth(0);
		chartPanel.setMinimumDrawHeight(0);
		chartPanel.setMaximumDrawWidth(Integer.MAX_VALUE);
		chartPanel.setMaximumDrawHeight(Integer.MAX_VALUE);
		
		return chartPanel;
	}

}
