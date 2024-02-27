package edu.washington.gs.maccoss.encyclopedia.utils.graphing;

import java.awt.Color;
import java.util.Optional;

public class CategoricalData {
	private final String category;
	private final float[] data;
	private final Optional<Color> color;
	
	public CategoricalData(String category, float[] data) {
		this.category = category;
		this.data = data;
		this.color = Optional.empty();
	}
	public CategoricalData(String category, float[] data, Color color) {
		this.category = category;
		this.data = data;
		this.color = Optional.of(color);
	}
	public String getCategory() {
		return category;
	}
	public float[] getData() {
		return data;
	}
	public Optional<Color> getColor() {
		return color;
	}
}
