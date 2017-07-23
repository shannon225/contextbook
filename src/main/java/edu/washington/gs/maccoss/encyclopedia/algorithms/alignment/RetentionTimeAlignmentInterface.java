package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;

public interface RetentionTimeAlignmentInterface {

	void plot(ArrayList<XYPoint> rts, Optional<File> saveFileSeed);

	float getYValue(float xrt);

	float getXValue(float yrt);

	float getProbabilityFitsModel(float actualRT, float modelRT);

	float getProbabilityFitsModel(float delta);

}