package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.font.TextAttribute;
import java.awt.geom.Rectangle2D;
import java.text.AttributedString;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.axis.LogAxis;
import org.jfree.chart.axis.LogTick;
import org.jfree.chart.axis.TickType;
import org.jfree.data.Range;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.TextAnchor;

import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.procedure.TDoubleObjectProcedure;

public class ExtendedLogAxis extends LogAxis {
	private static final long serialVersionUID = 1L;
	
	private final TDoubleObjectHashMap<String> automaticTickLabels=new TDoubleObjectHashMap<String>();

	public ExtendedLogAxis() {
		super();
	}

	public ExtendedLogAxis(String label) {
		super(label);
	}
    
    private double calculateValueNoINF(double log) {
        double result = calculateValue(log);
        if (Double.isInfinite(result)) {
            result = Double.MAX_VALUE;
        }
        if (result <= 0.0) {
            result = Double.MIN_VALUE;
        }
        return result;
    }

    /**
     * Returns a list of ticks for an axis at the top or bottom of the chart.
     *
     * @param g2  the graphics device ({@code null} not permitted).
     * @param dataArea  the data area ({@code null} not permitted).
     * @param edge  the edge ({@code null} not permitted).
     *
     * @return A list of ticks.
     */
    protected List refreshTicksHorizontal(Graphics2D g2, Rectangle2D dataArea,
            RectangleEdge edge) {

        Range range = getRange();
        List ticks = new ArrayList();
        Font tickLabelFont = getTickLabelFont();
        g2.setFont(tickLabelFont);
        TextAnchor textAnchor;
        if (edge == RectangleEdge.TOP) {
            textAnchor = TextAnchor.BOTTOM_CENTER;
        }
        else {
            textAnchor = TextAnchor.TOP_CENTER;
        }

        if (isAutoTickUnitSelection()) {
            selectAutoTickUnit(g2, dataArea, edge);
        }
        int minorTickCount = getTickUnit().getMinorTickCount();
        double unit = getTickUnit().getSize();
        double index = Math.ceil(calculateLog(getRange().getLowerBound()) 
                / unit);
        double start = index * unit;
        double end = calculateLog(getUpperBound());
        double current = start;
        boolean hasTicks = (getTickUnit().getSize() > 0.0)
                           && !Double.isInfinite(start);
        while (hasTicks && current <= end) {
            double v = calculateValueNoINF(current);
            if (range.contains(v)) {
                LogTick majorTick = getMajorTick(textAnchor, v);
                if (majorTick!=null) {
                	ticks.add(majorTick);
                }
            }
            // add minor ticks (for gridlines)
            double next = Math.pow(getBase(), current
                    + getTickUnit().getSize());
            for (int i = 1; i < minorTickCount; i++) {
                double minorV = v + i * ((next - v) / minorTickCount);
                if (range.contains(minorV)) {
                    LogTick minorTick = getMinorTick(textAnchor, minorV);
                    if (minorTick!=null) {
                    	ticks.add(minorTick);
                    }
                }
            }
            current = current + getTickUnit().getSize();
        }
        
        automaticTickLabels.forEachEntry(new TDoubleObjectProcedure<String>() {
        	public boolean execute(double v, String text) {
                if (range.contains(v)) {
                	AttributedString label = new AttributedString(text);
                	label.addAttributes(getTickLabelFont().getAttributes(), 0, text.length());
					ticks.add(new LogTick(TickType.MAJOR, v, label, textAnchor));
                }
        		return true;
        	};
		});
        return ticks;
    }

    /**
     * Returns a list of ticks for an axis at the left or right of the chart.
     *
     * @param g2  the graphics device ({@code null} not permitted).
     * @param dataArea  the data area ({@code null} not permitted).
     * @param edge  the edge that the axis is aligned to ({@code null} 
     *     not permitted).
     *
     * @return A list of ticks.
     */
    protected List refreshTicksVertical(Graphics2D g2, Rectangle2D dataArea,
            RectangleEdge edge) {

        Range range = getRange();
        List ticks = new ArrayList();
        Font tickLabelFont = getTickLabelFont();
        g2.setFont(tickLabelFont);
        TextAnchor textAnchor;
        if (edge == RectangleEdge.RIGHT) {
            textAnchor = TextAnchor.CENTER_LEFT;
        }
        else {
            textAnchor = TextAnchor.CENTER_RIGHT;
        }

        if (isAutoTickUnitSelection()) {
            selectAutoTickUnit(g2, dataArea, edge);
        }
        int minorTickCount = getTickUnit().getMinorTickCount();
        double unit = getTickUnit().getSize();
        double index = Math.ceil(calculateLog(getRange().getLowerBound()) 
                / unit);
        double start = index * unit;
        double end = calculateLog(getUpperBound());
        double current = start;
        boolean hasTicks = (getTickUnit().getSize() > 0.0)
                           && !Double.isInfinite(start);
        while (hasTicks && current <= end) {
            double v = calculateValueNoINF(current);
            if (range.contains(v)) {
                LogTick majorTick = getMajorTick(textAnchor, v);
                if (majorTick!=null) {
                	ticks.add(majorTick);
                }
            }
            // add minor ticks (for gridlines)
            double next = Math.pow(getBase(), current
                    + getTickUnit().getSize());
            for (int i = 1; i < minorTickCount; i++) {
                double minorV = v + i * ((next - v) / minorTickCount);
                if (range.contains(minorV)) {
                    LogTick minorTick = getMinorTick(textAnchor, minorV);
                    if (minorTick!=null) {
                    	ticks.add(minorTick);
                    }
                }
            }
            current = current + getTickUnit().getSize();
        }
        
        automaticTickLabels.forEachEntry(new TDoubleObjectProcedure<String>() {
        	public boolean execute(double current, String text) {
                double v = calculateValueNoINF(current);
                if (range.contains(v)) {
                	AttributedString label = new AttributedString(text);
                	label.addAttributes(getTickLabelFont().getAttributes(), 0, text.length());
					ticks.add(new LogTick(TickType.MAJOR, v, label, textAnchor));
                }
        		return true;
        	};
		});
        return ticks;
    }

	protected LogTick getMinorTick(TextAnchor textAnchor, double minorV) {
		return new LogTick(TickType.MINOR, minorV, null, textAnchor);
	}

	protected LogTick getMajorTick(TextAnchor textAnchor, double v) {
		return new LogTick(TickType.MAJOR, v, createTickLabel(v), textAnchor);
	}
	
	public void addAutomaticTicks(double v, String text) {
		automaticTickLabels.put(v, text);
	}
}
