package edu.washington.gs.maccoss.encyclopedia.utils.graphing;

import java.awt.Color;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;

public interface XYTraceInterface {

	Optional<Color> getColor();

	Optional<Float> getThickness();

	String getName();

	GraphType getType();

	Pair<double[], double[]> toArrays();

}