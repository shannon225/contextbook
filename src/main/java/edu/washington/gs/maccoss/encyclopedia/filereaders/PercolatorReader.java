package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.ArrayList;
import java.util.Map;

import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;

public class PercolatorReader {
	public static ArrayList<ScoredObject<String>> getPassingPeptides(File f, final float qValueThreshold) {
		final ArrayList<ScoredObject<String>> data=new ArrayList<ScoredObject<String>>();
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String psmID=row.get("PSMId");
				float qvalue=Float.parseFloat(row.get("q-value"));
				if (qvalue<qValueThreshold) {
					data.add(new ScoredObject<String>(qvalue, psmID));
				}
			}
		};

		TableParser.readTable(f, "\t", muscle);
		return data;
	}
}
