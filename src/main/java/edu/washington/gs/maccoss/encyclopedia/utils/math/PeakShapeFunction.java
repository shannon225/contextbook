package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;

public abstract class PeakShapeFunction implements Function {
	public static final int NUMBER_OF_STEPS=20;
	public static final double EPSILON=1e-6;

    // over-calculate trapezoids to integrate
    public double integrate(double left, double right) {
    	double width=right-left;
    	
    	double previousX=left;
    	double previousY=getYValue((float)previousX);
    	double integration=0.0f;
    	for (int i=1; i<=NUMBER_OF_STEPS; i++) {
    		double newX=left+i*width/NUMBER_OF_STEPS;
    		double newY=getYValue((float)newX);

    		double trap=(newX-previousX)*(newY+previousY)/2.0f;
			integration+=trap;
			
			previousX=newX;
			previousY=newY;
		}
    	return integration;
    }

	public XYPoint getApex(double left, double right) {
		double l=left, r=right;
		while ((r-l)>EPSILON) {
			double m1=l+(r-l)/3.0;
			double m2=r-(r-l)/3.0;
			double f1=getYValue((float)m1);
			double f2=getYValue((float)m2);
			if (f1<f2) l=m1;
			else r=m2;
		}
		
		double x=0.5*(l+r);

		// it probably won't, but keep a boundary if it wins
		double y=getYValue((float)x);
		double yL=getYValue((float)left);
		double yR=getYValue((float)right);
		if (yL>=y&&yL>=yR) return new XYPoint(left, yL);
		if (yR>=y&&yR>=yL) return new XYPoint(right, yR);
		return new XYPoint(x, y);
	}

	public ArrayList<XYPoint> interpolate() {
		ArrayList<XYPoint> points=getKnots();
		ArrayList<XYPoint> interpolated=new ArrayList<XYPoint>();
		
		XYPoint last=null;
		for (XYPoint xyPoint : points) {
			if (last==null) {
				last=xyPoint;
				interpolated.add(new XYPoint(xyPoint.x, getYValue((float)xyPoint.x)));
			} else {
				double delta=(xyPoint.x-last.x)/NUMBER_OF_STEPS;
				// interpolate the points in between
				for (int i=1; i<NUMBER_OF_STEPS; i++) {
					double newX=i*delta+last.x;
					double newY=getYValue((float)newX);
					interpolated.add(new XYPoint(newX, newY));
				}

				last=xyPoint;
				interpolated.add(new XYPoint(xyPoint.x, getYValue((float)xyPoint.x)));
			}
		}
		return interpolated;
	}
}
