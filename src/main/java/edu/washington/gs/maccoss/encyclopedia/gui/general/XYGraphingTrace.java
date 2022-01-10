package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.geom.Ellipse2D;
import java.util.Optional;

import org.jfree.chart.renderer.xy.AbstractXYItemRenderer;
import org.jfree.chart.renderer.xy.XYAreaRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;

public class XYGraphingTrace extends XYTrace {
	public XYGraphingTrace(XYTrace trace) {
		super(trace.getPoints(), trace.getType(), trace.getName(), trace.getColor(), trace.getThickness());
	}
	
	public XYGraphingTrace(XYTrace trace, GraphType type, String name) {
		super(trace.getPoints(), type, name);
	}
	
	public XYGraphingTrace(XYTrace trace, GraphType type, String name, Optional<Color> color, Optional<Float> thickness) {
		super(trace.getPoints(), type, name, color.orElse(null), thickness.orElse(null));
	}

	public XYSeries getSeries() {
		XYSeries series=new XYSeries(getName());
		for (XYPoint xy : getPoints()) {
			if (!Double.isNaN(xy.x)&&!Double.isNaN(xy.y)) {
				series.add(xy.x, xy.y);
			}
		}
		return series;
	}
	
	public AbstractXYItemRenderer getRenderer() {
		AbstractXYItemRenderer renderer;
		switch (getType()) {
		case area:
			renderer=new XYAreaRenderer();
			renderer.setSeriesStroke(0, new BasicStroke(2.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
			break;

		case line:
			renderer=new XYLineAndShapeRenderer();
			renderer.setSeriesStroke(0, new BasicStroke(getThickness().orElse(2.0f), BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
			((XYLineAndShapeRenderer) renderer).setBaseShapesVisible(false);

			break;

		case boldline:
			renderer=new XYLineAndShapeRenderer();
			renderer.setSeriesStroke(0, new BasicStroke(getThickness().orElse(5.0f), BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
			((XYLineAndShapeRenderer) renderer).setBaseShapesVisible(false);

			break;

		case squaredline:
			renderer=new XYLineAndShapeRenderer();
			renderer.setSeriesStroke(0, new BasicStroke(getThickness().orElse(5.0f), BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL));
			((XYLineAndShapeRenderer) renderer).setBaseShapesVisible(false);

			break;

		case dashedline:
			renderer=new XYLineAndShapeRenderer();
			Float thickness=getThickness().orElse(2.0f);
			if (thickness>5) {
				renderer.setSeriesStroke(0, new BasicStroke(thickness, BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 10.0f, new float[] {12.0f, 16.0f}, 0.0f));
			} else {
				renderer.setSeriesStroke(0, new BasicStroke(thickness, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 0.0f, new float[] {3.0f, 5.0f}, 0.0f));
			}
			((XYLineAndShapeRenderer) renderer).setDrawSeriesLineAsPath(true);
			((XYLineAndShapeRenderer) renderer).setBaseShapesVisible(false);

			break;

		case bigpoint:
			renderer=new XYLineAndShapeRenderer();
			renderer.setSeriesShape(0, new Ellipse2D.Double(-3, -3, 6, 6));
			((XYLineAndShapeRenderer) renderer).setBaseLinesVisible(false);

			break;

		case point:
			renderer=new XYLineAndShapeRenderer();
			renderer.setSeriesShape(0, new Ellipse2D.Double(-1.5, -1.5, 3, 3));
			((XYLineAndShapeRenderer) renderer).setBaseLinesVisible(false);

			break;

		case tinypoint:
			renderer=new XYLineAndShapeRenderer();
			renderer.setSeriesShape(0, new Ellipse2D.Double(-0.5, -0.5, 1, 1));
			((XYLineAndShapeRenderer) renderer).setBaseLinesVisible(false);

			break;

		case text:
			renderer=new XYLineAndShapeRenderer();
			renderer.setSeriesShape(0, new Ellipse2D.Double(-0.5, -0.5, 1, 1));
			((XYLineAndShapeRenderer) renderer).setBaseLinesVisible(false);

			break;

		case spectrum:
			renderer=new XYLineAndShapeRenderer();
			((XYLineAndShapeRenderer) renderer).setBaseShapesVisible(false);
			renderer.setBasePaint(Color.DARK_GRAY);

			break;

		default:
			throw new EncyclopediaException("unsupported graphing type!");
		}
		return renderer;
	}
}
