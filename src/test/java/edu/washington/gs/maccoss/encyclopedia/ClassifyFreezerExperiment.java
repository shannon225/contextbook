package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.util.ArrayList;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.math.randomforest.Record;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class ClassifyFreezerExperiment {
	public static ArrayList<Record> getRecordsFromTSV(File f) {
		final ArrayList<Record> data=new ArrayList<Record>();
		
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String sample=row.get("sample");
				float tic=Float.parseFloat(row.get("tic"));
				data.put(sample, tic);
			}
		};
		
		TableParser.parseTSV(f, muscle);

		return data;
	}
}
