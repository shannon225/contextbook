package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.Axis;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.entity.EntityCollection;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.CrosshairState;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotRenderingInfo;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.Outlier;
import org.jfree.chart.renderer.OutlierList;
import org.jfree.chart.renderer.OutlierListCollection;
import org.jfree.chart.renderer.category.BoxAndWhiskerRenderer;
import org.jfree.chart.renderer.category.CategoryItemRendererState;
import org.jfree.chart.renderer.xy.XYBoxAndWhiskerRenderer;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.statistics.BoxAndWhiskerCategoryDataset;
import org.jfree.data.statistics.BoxAndWhiskerItem;
import org.jfree.data.statistics.BoxAndWhiskerXYDataset;
import org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.text.TextUtilities;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.TextAnchor;

import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import gnu.trove.list.array.TFloatArrayList;

public class Boxplotter {
	private static Paint linePaint=Color.BLACK;
	private static Stroke darklineStroke=new BasicStroke(2.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL);
	private static Stroke lineStroke=new BasicStroke(0.5f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL);
	private static Stroke lightlineStroke=new BasicStroke(0.5f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL);
	private static Font lineFont=new Font("News Gothic MT", Font.PLAIN, 10);

	public static void main(String[] args) {
		XYBoxPlotterRenderer renderer=new XYBoxPlotterRenderer(100, false);

		final int seriesCount=3;
		final int categoryCount=4;
		final int entityCount=220;

		final NumberBoxAndWhiskerXYDataset dataset=new NumberBoxAndWhiskerXYDataset("Set");
		double min=Double.MAX_VALUE;
		double max=-Double.MAX_VALUE;
		int count=0;
		for (int i=0; i<seriesCount; i++) {
			for (int j=0; j<categoryCount; j++) {
				TFloatArrayList list=new TFloatArrayList();
				// add some values...
				for (int k=0; k<entityCount; k++) {
					list.add(7.0f+(float)Math.random()*6+i-j);
					list.add(9.0f+(float)Math.random()*2+i-j);
				}
				for (int k=0; k<entityCount/10; k++) {
					list.add(7.0f+(float)Math.random()*20+i-j);
				}
				for (float d : list.toArray()) {
					if (d<min) min=d;
					if (d>max) max=d;
				}
				count++;
				dataset.add(count, calculateBoxAndWhiskerStatistics(list.toArray()));
			}
		}

		final NumberAxis xAxis=new NumberAxis("Each Boxplot");
		final NumberAxis yAxis=new NumberAxis("Value");
		yAxis.setAutoRangeIncludesZero(false);
		final XYPlot plot=new XYPlot(dataset, xAxis, yAxis, renderer);

		final JFreeChart chart=new JFreeChart("Box-and-Whisker XY Demo", new Font("SansSerif", Font.BOLD, 14), plot, true);
		final ChartPanel chartPanel=new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(450, 270));

		Charter.launchChart(chartPanel, "Box-and-Whisker XY Demo");
	}
	
	public static void main2(String[] args) {
		CategoryBoxPlotterRenderer renderer=new CategoryBoxPlotterRenderer();

		final int seriesCount=3;
		final int categoryCount=4;
		final int entityCount=220;

		final DefaultBoxAndWhiskerCategoryDataset dataset=new DefaultBoxAndWhiskerCategoryDataset();
		double min=Double.MAX_VALUE;
		double max=-Double.MAX_VALUE;
		for (int i=0; i<seriesCount; i++) {
			for (int j=0; j<categoryCount; j++) {
				TFloatArrayList list=new TFloatArrayList();
				// add some values...
				for (int k=0; k<entityCount; k++) {
					list.add(7.0f+(float)Math.random()*6+i-j);
					list.add(9.0f+(float)Math.random()*2+i-j);
				}
				for (int k=0; k<entityCount/10; k++) {
					list.add(7.0f+(float)Math.random()*20+i-j);
				}
				for (float d : list.toArray()) {
					if (d<min) min=d;
					if (d>max) max=d;
				}
				dataset.add(calculateBoxAndWhiskerStatistics(list.toArray()), "Series "+i, " Type "+j);
			}
		}

		final CategoryAxis xAxis=new CategoryAxis("Type");
		final NumberAxis yAxis=new NumberAxis("Value");
		yAxis.setAutoRangeIncludesZero(false);
		final CategoryPlot plot=new CategoryPlot(dataset, xAxis, yAxis, renderer);

		final JFreeChart chart=new JFreeChart("Box-and-Whisker Demo", new Font("SansSerif", Font.BOLD, 14), plot, true);
		final ChartPanel chartPanel=new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(450, 270));

		Charter.launchChart(chartPanel, "Box-and-Whisker Demo");
	}

	public static BoxAndWhiskerItem calculateINFProtectedBoxAndWhiskerStatistics(float[] f, float positiveInf, float negativeInf, float nan) {
		return calculateINFProtectedBoxAndWhiskerStatistics(f, positiveInf, negativeInf, nan, false);
	}

	public static BoxAndWhiskerItem calculateINFProtectedBoxAndWhiskerStatistics(float[] f, float positiveInf, float negativeInf, float nan, boolean isLogarithmic) {
		float[] data=f.clone();

		for (int i = 0; i < data.length; i++) {
			if (Float.isNaN(data[i])) {
				data[i]=nan;
			} else if (Float.isInfinite(data[i])) {
				if (data[i]>0) {
					data[i]=positiveInf;
				} else {
					data[i]=negativeInf;
				}
			} else if (data[i]>positiveInf) {
				data[i]=positiveInf;
			} else if (data[i]<negativeInf) {
				data[i]=negativeInf;
			}
		}
		
		return calculateBoxAndWhiskerStatistics(data, isLogarithmic);
	}

	public static BoxAndWhiskerItem calculateBoxAndWhiskerStatistics(float[] f) {
		return calculateBoxAndWhiskerStatistics(f, false);
	}
	public static BoxAndWhiskerItem calculateBoxAndWhiskerStatistics(float[] f, boolean isLogarithmic) {
		float[] data=f.clone();
		if (isLogarithmic) {
			data=Log.log10(data);
		}
		float mean=General.mean(data);
		float median=QuickMedian.select(data, 0.5f);
		float q1=QuickMedian.select(data, 0.25f);
		float q3=QuickMedian.select(data, 0.75f);
		float lowerOutlierThreshold=QuickMedian.select(data, 0.05f);
		float upperOutlierThreshold=QuickMedian.select(data, 0.95f);

		double minRegularValue=Double.POSITIVE_INFINITY;
		double maxRegularValue=Double.NEGATIVE_INFINITY;
		double minOutlier=Double.POSITIVE_INFINITY;
		double maxOutlier=Double.NEGATIVE_INFINITY;
		List<Float> outliers=new ArrayList<Float>();
		for (float value : data) {
			if (value>upperOutlierThreshold) {
				if (isLogarithmic) {
					outliers.add((float)Math.pow(10, value));
				} else {
					outliers.add(value);
				}
				if (value>maxOutlier) {
					maxOutlier=value;
				}
			} else if (value<lowerOutlierThreshold) {
				if (isLogarithmic) {
					outliers.add((float)Math.pow(10, value));
				} else {
					outliers.add(value);
				}
				if (value<minOutlier) {
					minOutlier=value;
				}
			} else {
				minRegularValue=Math.min(minRegularValue, value);
				maxRegularValue=Math.max(maxRegularValue, value);
			}
			minOutlier=Math.min(minOutlier, minRegularValue);
			maxOutlier=Math.max(maxOutlier, maxRegularValue);
		}
		if (minRegularValue==Double.POSITIVE_INFINITY) {
			minRegularValue=General.min(data);
		}
		if (maxRegularValue==Double.NEGATIVE_INFINITY) {
			maxRegularValue=General.max(data);
		}
		if (minOutlier==Double.POSITIVE_INFINITY) {
			minOutlier=General.min(data);
		}
		if (maxOutlier==Double.NEGATIVE_INFINITY) {
			maxOutlier=General.max(data);
		}

		if (isLogarithmic) {
			return new BoxAndWhiskerItem(Double.valueOf(Math.pow(10, mean)), Double.valueOf(Math.pow(10, median)), Double.valueOf(Math.pow(10, q1)), Double.valueOf(Math.pow(10, q3)), Double.valueOf(Math.pow(10, minRegularValue)), 
					Double.valueOf(Math.pow(10, maxRegularValue)), Double.valueOf(Math.pow(10, minOutlier)), Double.valueOf(Math.pow(10, maxOutlier)), outliers);
		} else {
			return new BoxAndWhiskerItem(Double.valueOf(mean), Double.valueOf(median), Double.valueOf(q1), Double.valueOf(q3), Double.valueOf(minRegularValue), 
					Double.valueOf(maxRegularValue), Double.valueOf(minOutlier), Double.valueOf(maxOutlier), outliers);
		}
	}
	
	public static class XYBoxPlotterRenderer extends XYBoxAndWhiskerRenderer {
		private static final long serialVersionUID=1L;
		private final double barWidth;
		private final boolean isLogarithmic;
		
		public XYBoxPlotterRenderer(int barWidth, boolean isLogarithmic) {
			super();
			this.barWidth=barWidth;
			this.isLogarithmic=isLogarithmic;
		}

		public void drawVerticalItem(Graphics2D g2, Rectangle2D dataArea, PlotRenderingInfo info, XYPlot plot,
				ValueAxis domainAxis, ValueAxis rangeAxis, XYDataset dataset, int row, int column,
				CrosshairState crosshairState, int pass) {
	    	BoxAndWhiskerXYDataset bawDataset=(BoxAndWhiskerXYDataset)dataset;

	        Number x = bawDataset.getX(row, column);
	        double xx = domainAxis.valueToJava2D(x.doubleValue(), dataArea,
	                plot.getDomainAxisEdge())-barWidth/2.0;

			RectangleEdge rangeAxisEdge=plot.getRangeAxisEdge();
			Paint itemPaint=getItemPaint(row, column);
			Shape box=getVerticalItem(g2, dataArea, plot, domainAxis, rangeAxis, BoxAndWhiskerDataInterfaceGenerator.getInterface(bawDataset), row, column, xx, barWidth, rangeAxisEdge, itemPaint, isLogarithmic);
			
			EntityCollection entities = null;
	        if (info != null) {
	            entities = info.getOwner().getEntityCollection();
	            
	            if (box.intersects(dataArea)) {
	                Number yAverage = bawDataset.getMeanValue(row, column);
	                double yyAverage = rangeAxis.valueToJava2D(yAverage.doubleValue(), dataArea, rangeAxisEdge);
	                addEntity(entities, box, dataset, row, column, yyAverage, xx);
	            }
	        }
	    }
	}

	public static class CategoryBoxPlotterRenderer extends BoxAndWhiskerRenderer {
		private static final long serialVersionUID=1L;
		private final boolean isLogarithmic;

		public CategoryBoxPlotterRenderer() {
			this(false);
		}
		
		public CategoryBoxPlotterRenderer(boolean isLogarithmic) {
			super();
			this.isLogarithmic=isLogarithmic;
		}

		@Override
		public void drawVerticalItem(Graphics2D g2, CategoryItemRendererState state, Rectangle2D dataArea, CategoryPlot plot, CategoryAxis domainAxis, ValueAxis rangeAxis, CategoryDataset dataset,
				int row, int column) {

			BoxAndWhiskerCategoryDataset bawDataset=(BoxAndWhiskerCategoryDataset)dataset;

			double categoryEnd=domainAxis.getCategoryEnd(column, getColumnCount(), dataArea, plot.getDomainAxisEdge());
			double categoryStart=domainAxis.getCategoryStart(column, getColumnCount(), dataArea, plot.getDomainAxisEdge());
			double categoryWidth=categoryEnd-categoryStart;

			double xx=categoryStart;
			int seriesCount=getRowCount();
			int categoryCount=getColumnCount();

			if (seriesCount>1) {
				double seriesGap=dataArea.getWidth()*getItemMargin()/(categoryCount*(seriesCount-1));
				double usedWidth=(state.getBarWidth()*seriesCount)+(seriesGap*(seriesCount-1));
				// offset the start of the boxes if the total width used is
				// smaller than the category width
				double offset=(categoryWidth-usedWidth)/2;
				xx=xx+offset+(row*(state.getBarWidth()+seriesGap));
			} else {
				// offset the start of the box if the box width is smaller than
				// the
				// category width
				double offset=(categoryWidth-state.getBarWidth())/2;
				xx=xx+offset;
			}
			double barWidth=state.getBarWidth();
			RectangleEdge rangeAxisEdge=plot.getRangeAxisEdge();
			Paint itemPaint=getItemPaint(row, column);
			Shape box=getVerticalItem(g2, dataArea, plot, domainAxis, rangeAxis, BoxAndWhiskerDataInterfaceGenerator.getInterface(bawDataset), row, column, xx, barWidth, rangeAxisEdge, itemPaint, isLogarithmic);
			
			// collect entity and tool tip information...
			if (state.getInfo()!=null&&box!=null) {
				EntityCollection entities=state.getEntityCollection();
				if (entities!=null) {
					addItemEntity(entities, dataset, row, column, box);
				}
			}
		}

		
	}
	public static Shape getVerticalItem(Graphics2D g2, Rectangle2D dataArea, Plot plot, Axis domainAxis, ValueAxis rangeAxis, BoxAndWhiskerDataInterface bawDataset,
			int row, int column, double xx, double barWidth, RectangleEdge rangeAxisEdge, Paint itemPaint, boolean isLogarithmic) {

		double yyOutlier;
		Paint localLinePaint=linePaint;
		if (itemPaint instanceof Color) {
			//localLinePaint=((Color)itemPaint).darker();
		}
		g2.setPaint(localLinePaint);
		g2.setStroke(lineStroke);
		g2.setFont(lineFont);

		double aRadius=0; // average radius

		Number yMedian=bawDataset.getMedianValue(row, column);
		Number yQ1=bawDataset.getQ1Value(row, column);
		Number yQ3=bawDataset.getQ3Value(row, column);
		Number yMax=bawDataset.getMaxRegularValue(row, column);
		Number yMin=bawDataset.getMinRegularValue(row, column);
		List yOutliers=bawDataset.getOutliers(row, column);
		
		Shape box=null;
		if (yQ1!=null&&yQ3!=null&&yMax!=null&&yMin!=null) {

			double yyQ1=rangeAxis.valueToJava2D(yQ1.doubleValue(), dataArea, rangeAxisEdge);
			double yyQ3=rangeAxis.valueToJava2D(yQ3.doubleValue(), dataArea, rangeAxisEdge);
			double yyMax=rangeAxis.valueToJava2D(yMax.doubleValue(), dataArea, rangeAxisEdge);
			double yyMin=rangeAxis.valueToJava2D(yMin.doubleValue(), dataArea, rangeAxisEdge);
			double xxmid=xx+barWidth/2.0;

			// draw the upper shadow...
			g2.draw(new Line2D.Double(xxmid, yyMax, xxmid, yyQ3));
			g2.draw(new Line2D.Double(xx+barWidth/4, yyMax, xx+barWidth*3/4, yyMax));

			// draw the lower shadow...
			g2.draw(new Line2D.Double(xxmid, yyMin, xxmid, yyQ1));
			g2.draw(new Line2D.Double(xx+barWidth/4, yyMin, xx+barWidth*3/4, yyMin));

			// draw the body...
			g2.setPaint(itemPaint);
			box=new Rectangle2D.Double(xx, Math.min(yyQ1, yyQ3), barWidth, Math.abs(yyQ1-yyQ3));
			g2.fill(box);
			g2.setStroke(lineStroke);
			g2.setPaint(localLinePaint);
			g2.draw(box);
		}

		g2.setPaint(localLinePaint);

		// draw median...
		if (yMedian!=null) {
			g2.setStroke(darklineStroke);
			double yyMedian=rangeAxis.valueToJava2D(yMedian.doubleValue(), dataArea, rangeAxisEdge);
			g2.draw(new Line2D.Double(xx, yyMedian, xx+barWidth, yyMedian));
		}

		// draw yOutliers...
		double maxAxisValue=rangeAxis.valueToJava2D(rangeAxis.getUpperBound(), dataArea, rangeAxisEdge)+aRadius;
		double minAxisValue=rangeAxis.valueToJava2D(rangeAxis.getLowerBound(), dataArea, rangeAxisEdge)-aRadius;

		g2.setStroke(lightlineStroke);

		// draw outliers
		double oRadius=barWidth/3; // outlier radius
		List outliers=new ArrayList();
		OutlierListCollection outlierListCollection=new OutlierListCollection();

		// From outlier array sort out which are outliers and put these into
		// a
		// list If there are any farouts, set the flag on the
		// OutlierListCollection
		double upperLogCorrected=isLogarithmic?Log.log10(rangeAxis.getUpperBound()):rangeAxis.getUpperBound();
		double lowerLogCorrected=isLogarithmic?Log.log10(rangeAxis.getLowerBound()):rangeAxis.getLowerBound();
		double fontHeight=0;//(g2.getFontMetrics().getAscent()/(maxAxisValue-minAxisValue))*(upperLogCorrected-lowerLogCorrected);
		double minOutlier=lowerLogCorrected-fontHeight;
		double maxOutlier=upperLogCorrected+fontHeight;
		double minRegular=yMin.doubleValue();
		double maxRegular=yMax.doubleValue();
		
		if (yOutliers!=null) {
			int aboveOutliers=0;
			int belowOutliers=0;
			
			for (int i=0; i<yOutliers.size(); i++) {
				Number number=(Number)yOutliers.get(i);
				double outlier=number.doubleValue();

				if (outlier>maxOutlier) {
					outlierListCollection.setHighFarOut(true);
					aboveOutliers++;
				} else if (outlier<minOutlier) {
					outlierListCollection.setLowFarOut(true);
					belowOutliers++;
				} else if (outlier>maxRegular) {
					yyOutlier=rangeAxis.valueToJava2D(outlier, dataArea, rangeAxisEdge);
					float random=0.5f;//RandomGenerator.random(number.hashCode());
					outliers.add(new Outlier(xx+barWidth/4.0+random*barWidth/2.0, yyOutlier, oRadius));
				} else if (outlier<minRegular) {
					yyOutlier=rangeAxis.valueToJava2D(outlier, dataArea, rangeAxisEdge);
					float random=0.5f;//RandomGenerator.random(number.hashCode());
					outliers.add(new Outlier(xx+barWidth/4.0+random*barWidth/2.0, yyOutlier, oRadius));
				}
			}

			// Process outliers. Each outlier is either added to the
			// appropriate outlier list or a new outlier list is made
			for (Object outlier : outliers) {
				outlierListCollection.add((Outlier)outlier);
			}

			for (Iterator iterator=outlierListCollection.iterator(); iterator.hasNext();) {
				OutlierList list=(OutlierList)iterator.next();
				Outlier outlier=list.getAveragedOutlier();
				Point2D point=outlier.getPoint();

				Ellipse2D dot=new Ellipse2D.Double(point.getX()+oRadius/2.0, point.getY()+oRadius/2.0, oRadius, oRadius);
				g2.draw(dot);
			}
			// draw farout indicators
			if (outlierListCollection.isHighFarOut()) {
				TextUtilities.drawAlignedString("+"+aboveOutliers, g2, (float)(xx+barWidth/2.0), 
						(float)maxAxisValue, TextAnchor.TOP_CENTER);
			}

			if (outlierListCollection.isLowFarOut()) {
				TextUtilities.drawAlignedString("+"+belowOutliers, g2, (float)(xx+barWidth/2.0), 
						(float)minAxisValue, TextAnchor.BOTTOM_CENTER);
			}
		}
		
		return box;
	}
}
