package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.awt.Color;

public interface Ion {

	double getMass();

	String getName();

	Color getColor();
	
	byte getIndex();

	IonType getType();

}