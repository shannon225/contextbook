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
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.entity.EntityCollection;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.renderer.Outlier;
import org.jfree.chart.renderer.OutlierList;
import org.jfree.chart.renderer.OutlierListCollection;
import org.jfree.chart.renderer.category.BoxAndWhiskerRenderer;
import org.jfree.chart.renderer.category.CategoryItemRendererState;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.statistics.BoxAndWhiskerCategoryDataset;
import org.jfree.data.statistics.BoxAndWhiskerItem;
import org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;
import org.jfree.text.TextUtilities;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.TextAnchor;

import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import gnu.trove.list.array.TFloatArrayList;

public class Boxplotter {
	private static Paint linePaint=Color.DARK_GRAY;
	private static Stroke lineStroke=new BasicStroke(2.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
	private static Stroke lightlineStroke=new BasicStroke(1.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
	private static Font lineFont=new Font("News Gothic MT", Font.PLAIN, 10);

	public static void main(String[] args) {
		BoxPlotterRenderer renderer=new BoxPlotterRenderer();

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
					list.add(7.0f+(float)Math.random()*6);
					list.add(9.0f+(float)Math.random()*2);
				}
				for (int k=0; k<entityCount/10; k++) {
					list.add(7.0f+(float)Math.random()*20);
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

	public static BoxAndWhiskerItem calculateBoxAndWhiskerStatistics(float[] f) {
		float[] data=f.clone();
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
		for (float value : f) {
			if (value>upperOutlierThreshold) {
				outliers.add(value);
				if (value>maxOutlier) {
					maxOutlier=value;
				}
			} else if (value<lowerOutlierThreshold) {
				outliers.add(value);
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

		return new BoxAndWhiskerItem(new Double(mean), new Double(median), new Double(q1), new Double(q3), new Double(minRegularValue), new Double(maxRegularValue), new Double(minOutlier),
				new Double(maxOutlier), outliers);
	}

	public static class BoxPlotterRenderer extends BoxAndWhiskerRenderer {
		private static final long serialVersionUID=1L;

		public BoxPlotterRenderer() {
			super();
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

			double yyOutlier;
			Paint itemPaint=getItemPaint(row, column);
			g2.setPaint(linePaint);
			g2.setStroke(lineStroke);
			g2.setFont(lineFont);

			double aRadius=0; // average radius

			RectangleEdge location=plot.getRangeAxisEdge();

			Number yQ1=bawDataset.getQ1Value(row, column);
			Number yQ3=bawDataset.getQ3Value(row, column);
			Number yMax=bawDataset.getMaxRegularValue(row, column);
			Number yMin=bawDataset.getMinRegularValue(row, column);
			Shape box=null;
			if (yQ1!=null&&yQ3!=null&&yMax!=null&&yMin!=null) {

				double yyQ1=rangeAxis.valueToJava2D(yQ1.doubleValue(), dataArea, location);
				double yyQ3=rangeAxis.valueToJava2D(yQ3.doubleValue(), dataArea, location);
				double yyMax=rangeAxis.valueToJava2D(yMax.doubleValue(), dataArea, location);
				double yyMin=rangeAxis.valueToJava2D(yMin.doubleValue(), dataArea, location);
				double xxmid=xx+state.getBarWidth()/2.0;

				// draw the upper shadow...
				g2.draw(new Line2D.Double(xxmid, yyMax, xxmid, yyQ3));
				g2.draw(new Line2D.Double(xx+state.getBarWidth()/4, yyMax, xx+state.getBarWidth()*3/4, yyMax));

				// draw the lower shadow...
				g2.draw(new Line2D.Double(xxmid, yyMin, xxmid, yyQ1));
				g2.draw(new Line2D.Double(xx+state.getBarWidth()/4, yyMin, xx+state.getBarWidth()*3/4, yyMin));

				// draw the body...
				g2.setPaint(itemPaint);
				box=new Rectangle2D.Double(xx, Math.min(yyQ1, yyQ3), state.getBarWidth(), Math.abs(yyQ1-yyQ3));
				g2.fill(box);
				g2.setStroke(lineStroke);
				g2.setPaint(linePaint);
				g2.draw(box);
			}

			g2.setPaint(linePaint);

			// draw median...
			Number yMedian=bawDataset.getMedianValue(row, column);
			if (yMedian!=null) {
				double yyMedian=rangeAxis.valueToJava2D(yMedian.doubleValue(), dataArea, location);
				g2.draw(new Line2D.Double(xx, yyMedian, xx+state.getBarWidth(), yyMedian));
			}

			// draw yOutliers...
			double maxAxisValue=rangeAxis.valueToJava2D(rangeAxis.getUpperBound(), dataArea, location)+aRadius;
			double minAxisValue=rangeAxis.valueToJava2D(rangeAxis.getLowerBound(), dataArea, location)-aRadius;

			g2.setStroke(lightlineStroke);

			// draw outliers
			double oRadius=state.getBarWidth()/4; // outlier radius
			List outliers=new ArrayList();
			OutlierListCollection outlierListCollection=new OutlierListCollection();

			// From outlier array sort out which are outliers and put these into
			// a
			// list If there are any farouts, set the flag on the
			// OutlierListCollection
			List yOutliers=bawDataset.getOutliers(row, column);
			double fontHeight=(g2.getFontMetrics().getAscent()/(maxAxisValue-minAxisValue))*(rangeAxis.getUpperBound()-rangeAxis.getLowerBound());
			double minOutlier=rangeAxis.getLowerBound()-fontHeight;
			double maxOutlier=rangeAxis.getUpperBound()+fontHeight;
			double minRegular=bawDataset.getMinRegularValue(row, column).doubleValue();
			double maxRegular=bawDataset.getMaxRegularValue(row, column).doubleValue();
			
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
						yyOutlier=rangeAxis.valueToJava2D(outlier, dataArea, location);
						float random=RandomGenerator.random(number.hashCode());
						outliers.add(new Outlier(xx+state.getBarWidth()/4.0+random*state.getBarWidth()/2.0, yyOutlier, oRadius));
					} else if (outlier<minRegular) {
						yyOutlier=rangeAxis.valueToJava2D(outlier, dataArea, location);
						float random=RandomGenerator.random(number.hashCode());
						outliers.add(new Outlier(xx+state.getBarWidth()/4.0+random*state.getBarWidth()/2.0, yyOutlier, oRadius));
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

					Ellipse2D dot=new Ellipse2D.Double(point.getX()+oRadius/2, point.getY(), oRadius, oRadius);
					g2.draw(dot);
				}
				// draw farout indicators
				if (outlierListCollection.isHighFarOut()) {
					TextUtilities.drawAlignedString("+"+aboveOutliers, g2, (float)(xx+state.getBarWidth()/2.0), 
							(float)maxAxisValue, TextAnchor.TOP_CENTER);
				}

				if (outlierListCollection.isLowFarOut()) {
					TextUtilities.drawAlignedString("+"+belowOutliers, g2, (float)(xx+state.getBarWidth()/2.0), 
							(float)minAxisValue, TextAnchor.BOTTOM_CENTER);
				}
			}
			// collect entity and tool tip information...
			if (state.getInfo()!=null&&box!=null) {
				EntityCollection entities=state.getEntityCollection();
				if (entities!=null) {
					addItemEntity(entities, dataset, row, column, box);
				}
			}
		}
	}
}
