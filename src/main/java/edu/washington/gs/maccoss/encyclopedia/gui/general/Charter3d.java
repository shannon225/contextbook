package edu.washington.gs.maccoss.encyclopedia.gui.general;

//import org.jzy3d.analysis.AbstractAnalysis;
//import org.jzy3d.analysis.AnalysisLauncher;
//import org.jzy3d.chart.factories.AWTChartComponentFactory;
//import org.jzy3d.colors.Color;
//import org.jzy3d.colors.ColorMapper;
//import org.jzy3d.colors.colormaps.ColorMapRainbow;
//import org.jzy3d.maths.Range;
//import org.jzy3d.plot3d.builder.Builder;
//import org.jzy3d.plot3d.builder.Mapper;
//import org.jzy3d.plot3d.builder.concrete.OrthonormalGrid;
//import org.jzy3d.plot3d.primitives.Shape;
//import org.jzy3d.plot3d.rendering.canvas.Quality;
//
//import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
//
//public class Charter3d extends AbstractAnalysis {
//	private final Mapper mapper;
//	private final Range xrange;
//	private final Range yrange;
//	private final int steps;
//	
//	public static void main(String[] args) {
//        // Define a function to plot
//        Mapper mapper = new Mapper() {
//            @Override
//            public double f(double x, double y) {
//                return x * Math.sin(x * y);
//            }
//        };
//        plot(mapper, new Range(-3, 3), new Range(-3, 3), 80);
//	}
//	
//	public static Color getColor(java.awt.Color c) {
//		return new Color(c.getRed(), c.getGreen(), c.getBlue(), c.getAlpha());
//	}
//	
//	private static final float FACTOR=0.7f;
//	public static Color getDarker(Color c) {
//        return new Color(Math.max(c.r*FACTOR, 0.0f),
//                         Math.max(c.g*FACTOR, 0.0f),
//                         Math.max(c.b*FACTOR, 0.0f),
//                         c.a);
//	}
//
//	public Charter3d(Mapper mapper, Range xrange, Range yrange, int steps) {
//		this.mapper=mapper;
//		this.xrange=xrange;
//		this.yrange=yrange;
//		this.steps=steps;
//	}
//
//
//	public static void plot(Mapper mapper, Range xrange, Range yrange, int steps) {
//		try {
//			AnalysisLauncher.open(new Charter3d(mapper, xrange, yrange, steps));
//		} catch (Exception e) {
//			throw new EncyclopediaException("Unexpected plotting error", e);
//		}
//    }
//
//    @Override
//    public void init() {
//        // Create the object to represent the function over the given range.
//        final Shape surface = Builder.buildOrthonormal(new OrthonormalGrid(xrange, steps, yrange, steps), mapper);
//        surface.setColorMapper(new ColorMapper(new ColorMapRainbow(), surface.getBounds().getZmin(), surface.getBounds().getZmax(), new Color(1, 1, 1, 1f)));
//        surface.setFaceDisplayed(true);
//        
//        surface.setWireframeDisplayed(false);
//        surface.setWireframeWidth(0.5f);
//        surface.setWireframeColor(Color.BLACK);
//
//        // Create a chart
//        chart = AWTChartComponentFactory.chart(Quality.Advanced, getCanvasType());
//        chart.getScene().getGraph().add(surface);
//    }
//}
