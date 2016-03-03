package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.util.ArrayList;
import java.util.Map;

import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.math.randomforest.Record;

public class ClassifyFreezerExperiment {
	public static ArrayList<Record> getRecordsFromTSV(File f) {
		final ArrayList<Record> data=new ArrayList<Record>();
		
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String sample=row.get("sample");
				float tic=Float.parseFloat(row.get("tic"));
				// FIXME data.put(sample, tic);
			}
		};
		
		TableParser.parseTSV(f, muscle);

		return data;
	}
}
