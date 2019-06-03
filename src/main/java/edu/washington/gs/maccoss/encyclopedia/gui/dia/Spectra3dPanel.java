package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.awt.BorderLayout;
import java.awt.Component;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.JPanel;

import org.jzy3d.chart.Chart;
import org.jzy3d.chart.factories.AWTChartComponentFactory;
import org.jzy3d.chart.factories.IChartComponentFactory.Toolkit;
import org.jzy3d.colors.Color;
import org.jzy3d.maths.Coord3d;
import org.jzy3d.plot3d.primitives.AbstractDrawable;
import org.jzy3d.plot3d.primitives.LineStrip;
import org.jzy3d.plot3d.primitives.Point;
import org.jzy3d.plot3d.primitives.Polygon;
import org.jzy3d.plot3d.primitives.axes.layout.renderers.ITickRenderer;
import org.jzy3d.plot3d.rendering.canvas.Quality;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter3d;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;

public class Spectra3dPanel extends JPanel {
	private static final long serialVersionUID=1L;
	private static final Color[] colors=new Color[] {new Color(102, 0, 153), new Color(204, 0, 0), new Color(0, 0, 255), new Color(0, 153, 0), new Color(255, 102, 0), new Color(102, 51, 0),
			new Color(152, 50, 203), new Color(255, 102, 102), new Color(51, 204, 255), new Color(0, 255, 51), new Color(255, 255, 0), new Color(255, 153, 0), new Color(153, 102, 0),};

	private final ArrayList<AbstractDrawable> peaks;
	private final Chart chart;

	public Spectra3dPanel(ArrayList<? extends Spectrum> entries, FragmentIon[] ions, MassTolerance tolerance) {
		super(new BorderLayout());

		float maxIntensity=0.0f;
		for (Spectrum s : entries) {
			maxIntensity=Math.max(maxIntensity, General.max(s.getIntensityArray()));
		}
		float multiplier=1.0f;
		if (maxIntensity>1000) {
			multiplier=(float)Math.pow(10, Math.floor(Log.log10(maxIntensity)));
		}
		
		System.out.println("test1");
		this.peaks=generatePeaks(entries, ions, tolerance, multiplier);
		System.out.println("test2");
		this.chart=generateChart(multiplier);
		System.out.println("test3: "+peaks.size()+" peaks");
		
		for (AbstractDrawable peak : peaks) {
			this.chart.getScene().add(peak);	
		}
		System.out.println("test4");

//		AWTCameraMouseController mouse=new AWTCameraMouseController();
//		chart.addController(mouse);
//
//		mouse.addControllerEventListener(new ControllerEventListener() {
//			public void controllerEventFired(ControllerEvent e) {
//				if (e.getType()==ControllerType.ROTATE) {
//					System.out.println("Mouse[VIEWPOINT]:"+chart.getViewPoint());
//				}
//			}
//		});

		this.add((Component)chart.getCanvas(), BorderLayout.CENTER);
	}
	
	public Chart getChart() {
		return chart;
	}

	private static final DecimalFormat xdf=new DecimalFormat("#");
	private static final DecimalFormat zdf=new DecimalFormat("0.0");
	public static Chart generateChart(float multiplier) {
		AWTChartComponentFactory factory=new AWTChartComponentFactory();
		Chart chart=factory.newChart(Quality.Nicest, Toolkit.awt.name());
		chart.getView().setMaximized(true);

		chart.getAxeLayout().setXAxeLabel("M/Z");
		chart.getAxeLayout().setYAxeLabel("Retention Time (min)");
		if (multiplier>1) {
			int order=Math.round(Log.log10(multiplier));
			chart.getAxeLayout().setZAxeLabel("Intensity (10e"+order+")");
		} else {
			chart.getAxeLayout().setZAxeLabel("Intensity");
		}
		ITickRenderer xrend = new ITickRenderer() {
			@Override
			public String format(double arg0) {
		        return xdf.format(arg0);
			}
		};
		ITickRenderer zrend = new ITickRenderer() {
			@Override
			public String format(double arg0) {
		        return zdf.format(arg0);
			}
		};
		chart.getAxeLayout().setXTickRenderer(xrend);
		chart.getAxeLayout().setZTickRenderer(zrend);
		chart.getAxeLayout().setYTickRenderer(zrend);
		chart.getAxeLayout().setGridColor(new Color(0, 0, 0, 0.5f));
		chart.getView().setViewPoint(new Coord3d(-1.3, 0.5, 0));
		
		return chart;
	}

	public static ArrayList<AbstractDrawable> generatePeaks(ArrayList<? extends Spectrum> entries, FragmentIon[] ions, MassTolerance tolerance, float multiplier) {
		ArrayList<AbstractDrawable> peaks=new ArrayList<>();
		HashMap<FragmentIon, Coord3d> lastPoints=new HashMap<>();
		
		for (Spectrum s : entries) {
			HashMap<FragmentIon, Coord3d> newPoints=new HashMap<>();
			double[] masses=s.getMassArray();
			float[] intensities=s.getIntensityArray();
			float rtInMin=s.getScanStartTime()/60f;
			
			for (int i=0; i<masses.length; i++) {
				Color color=Color.GRAY;
				FragmentIon ion=null;
				for (int j=0; j<ions.length; j++) {
					if (tolerance.equals(ions[j].mass, masses[i])) {
						//color=colors[j%colors.length];
						color=Charter3d.getColor(RandomGenerator.randomColor(ions[j].toString().hashCode()));
						if (((int)ions[j].mass)==150||((int)ions[j].mass)==133) color=Color.BLUE;
						ion=ions[j];
					}
				}
				Coord3d bottom=new Coord3d(masses[i], rtInMin, 0.0);
				Coord3d top=new Coord3d(masses[i], rtInMin, intensities[i]/multiplier);

				if (ion==null) {
					LineStrip peak=new LineStrip(new Point(bottom), new Point(top));
					peak.setWireframeColor(color);
					peaks.add(peak);
				}
				
				if (ion!=null) {
					newPoints.put(ion, top);
					Coord3d lastTop=lastPoints.get(ion);
					if (lastTop!=null) {
						Polygon cell=new Polygon();
						cell.add(new Point(new Coord3d(lastTop.x, lastTop.y, 0.0)));
						cell.add(new Point(lastTop));
						cell.add(new Point(top));
						cell.add(new Point(bottom));
						cell.setWireframeColor(Charter3d.getDarker(color));
						cell.setColor(color.alphaSelf(0.5f));
						cell.setFaceDisplayed(true);
						cell.setWireframeDisplayed(true);
						
						peaks.add(cell);
						
					}
				}
			}
			lastPoints=newPoints;
		}
		return peaks;
	}
}
