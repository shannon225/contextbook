package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearInterpolatedFunction;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedianDouble;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RTProbabilityModel;
import gnu.trove.list.array.TDoubleArrayList;

public class MedianRetentionTimeFilter extends AbstractRetentionTimeFilter {
	public static MedianRetentionTimeFilter getFilter(ArrayList<XYPoint> rts) {
		return getFilter(rts, RT_STRING, "Retention Time (min)");
	}
	
	public static MedianRetentionTimeFilter getFilter(ArrayList<XYPoint> rts, String xAxis, String yAxis) {
		Collections.sort(rts);
		int binsize=Math.round(Math.max(rts.size()/40.0f, 20));
		int delta=Math.round(rts.size()/(Math.max(rts.size()/40.0f, 20)*2-1));
		
		ArrayList<XYPoint> smoothedRTs=new ArrayList<>();
		final int n=rts.size()-binsize;
		for (int i=0; i<n; i+=delta) {
			TDoubleArrayList xs=new TDoubleArrayList();
			TDoubleArrayList ys=new TDoubleArrayList();
			for (int j=0; j<binsize; j++) {
				xs.add(rts.get(i+j).x);
				ys.add(rts.get(i+j).y);
			}
			final XYPoint xy=getMedianReference(xs.toArray(), ys.toArray());
			//System.out.println(i+"\t"+xy+"\t"+binsize);
			smoothedRTs.add(xy);
		}
		System.out.println(smoothedRTs.size()+" number of bins from N="+rts.size()+"!");
		
		Function rtWarper=new LinearInterpolatedFunction(smoothedRTs);
		
		Optional<RTProbabilityModel> model=Optional.of(RetentionTimeFilter.generateMixtureModel(rts, rtWarper));
		return new MedianRetentionTimeFilter(rtWarper, model, xAxis, yAxis);
	}

	private MedianRetentionTimeFilter(Function rtWarper, Optional<RTProbabilityModel> model, String xAxis, String yAxis) {
		super(rtWarper, model, xAxis, yAxis);
	}

	public static XYPoint getMedianReference(double[] v1, double[] v2) {
		double median1=QuickMedianDouble.median(v1);
		double median2=QuickMedianDouble.median(v2);

		return new XYPoint(median1, median2);
	}
}
