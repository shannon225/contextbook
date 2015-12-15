package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.Collections;
import java.util.List;
import java.util.Random;

public class RandomGenerator {
	/**
	 * rand() from ANSI C, should be faster (and easier to control) than java.util.Random, which uses Longs
	 * @param seed (Please don't use 0 as a seed!)
	 * @return
	 */
	public static int randomInt(int seed) {
		seed=seed*1103515245+12345;
		return seed%2147483647;
	}
	public static float random(int seed) {
		return floatFromRandomInt(randomInt(seed));
	}
	public static float floatFromRandomInt(int random) {
		return ((random/(float)2147483647)+1f)/2f;
	}
	public static float randomIndex(int length, int seed) {
		return (int)(length*floatFromRandomInt(randomInt(seed)));
	}

	public static void shuffle(List<?> list, int seed) {
		Collections.shuffle(list, new Random(seed));
    }
}
