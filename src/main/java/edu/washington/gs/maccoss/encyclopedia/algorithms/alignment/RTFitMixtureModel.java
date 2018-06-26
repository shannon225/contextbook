package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearInterpolatedFunction;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RTProbabilityModel;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Gaussian;
import gnu.trove.list.array.TFloatArrayList;

public class RTFitMixtureModel implements RTProbabilityModel {
	private final Function stdevFunction;
	
	public RTFitMixtureModel(ArrayList<XYPoint> rts, Function warper) {
		Collections.sort(rts);
		int binsize=(int)Math.round(rts.size()/Log.log2((float)rts.size())); //(int)(Math.sqrt(rts.size())); 
		TFloatArrayList deltas=new TFloatArrayList();
		for (int i=0; i<rts.size(); i++) {
			XYPoint xyPoint=rts.get(i);
			float delta=getDelta(warper, (float)xyPoint.y, (float)xyPoint.x);
			deltas.add(delta);
		}
		
		int half=binsize/2;
		
		ArrayList<XYPoint> stdevs=new ArrayList<>();
		for (int i=0; i<rts.size(); i++) {
			int start=Math.max(0, i-half);
			int stop=start+binsize;
			if (stop>=rts.size()) {
				stop=rts.size()-1;
				start=Math.max(0, stop-binsize);
			}
			
			TFloatArrayList localDeltas=new TFloatArrayList();
			for (int j=start; j<stop; j++) {
				localDeltas.add(deltas.get(j));
			}
			float[] data=localDeltas.toArray();
			float outside2stdevs=1.0f-0.6826f; // estimates based on 2x stdevs of data (middle 2/3rds)
			float stdev=(QuickMedian.select(data, 1.0f-outside2stdevs/2.0f)-QuickMedian.select(data, outside2stdevs/2.0f))/2.0f; //QuickMedian.iqr(localDeltas.toArray())/1.349f;
			
			
			stdevs.add(new XYPoint(rts.get(i).x, stdev));
		}
		
		stdevFunction=new LinearInterpolatedFunction(stdevs);
	}

	@Override
	public float getProbability(float retentionTime, float delta) {
		float stdev=stdevFunction.getYValue(retentionTime);
		return (float)(new Gaussian(0, stdev, 1.0f).getPDF(delta));
	}


	private static float getDelta(Function rtWarper, float actualRT, float modelRT) {
		float one=actualRT-rtWarper.getYValue(modelRT);
		float two=rtWarper.getXValue(actualRT)-modelRT;
		if (Math.abs(one)<Math.abs(two)) {
			return one;
		} else {
			return two;
		}
	}
}
