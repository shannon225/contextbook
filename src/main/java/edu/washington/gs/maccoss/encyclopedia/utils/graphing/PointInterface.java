package edu.washington.gs.maccoss.encyclopedia.utils.graphing;

public interface PointInterface extends Comparable<PointInterface> {

	public abstract double getX();

	public abstract double getY();

	public abstract int compareTo(PointInterface o);

}