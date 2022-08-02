package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;

public class RTReplacerForCE {
	public static void main(String[] args) throws Exception {
		File reference=new File("/Users/searle.30/Downloads/pan_human_library.dlib");
		File save=new File("/Users/searle.30/Downloads/pan_human_library_CE.dlib");

		LibraryFile referenceLibrary=new LibraryFile();
		referenceLibrary.openFile(reference);
		
		
		ArrayList<LibraryEntry> saveEntries=new ArrayList<>();
		for (LibraryEntry entry : referenceLibrary.getAllEntries(false, new AminoAcidConstants())) {
			byte Q=getQ(entry.getPeptideSeq());
			double M=MassConstants.getPeptideMass(entry.getPrecursorMZ(), Q);
			
			double ce = getCE(Q, M);
			saveEntries.add(entry.updateRetentionTime(60*(float)ce));
		}

		LibraryFile saveLibrary=new LibraryFile();
		saveLibrary.openFile();
		saveLibrary.addEntries(saveEntries);
		saveLibrary.saveAsFile(save);

	}

	private static byte getQ(String seq) {
		byte q=1;
		for(char c : seq.toCharArray()) {
			if ('K'==c) q++;
			if ('H'==c) q++;
			if ('R'==c) q++;
		}
		return q;
	}

	private static double getCE(float Q, double M) {
		double ce=Math.pow(M, 0.411)/Math.log(1+0.35*Q);
		return ce/5-2;
	}
	
	public static void main2(String[] args) throws Exception {
		File input=new File("/Users/searle.30/Downloads/RTFitFirst_CE/ce_peptides.txt");
		

		XYTraceInterface t1=new XYTrace(getData(input, 1, 0, 0), GraphType.tinypoint, "peptides", new Color(0f, 0f, 1f, 0.1f), 1.0f);
		//XYTraceInterface t2=new XYTrace(getData(input, 2, 0, 0), GraphType.tinypoint, "peptides", new Color(0f, 0.5f, 1f), 1.0f);

		//XYTraceInterface t2=new XYTrace(getData(input, 0, 1, 0), GraphType.tinypoint, "peptides", new Color(1f, 0f, 0f), 1.0f);
		//XYTraceInterface t2=new XYTrace(getData(input, 0, 2, 0), GraphType.tinypoint, "peptides", new Color(1f, 0f, 0f), 1.0f);
		
		//XYTraceInterface t2=new XYTrace(getData(input, 0, 1, 1), GraphType.tinypoint, "peptides", new Color(1f, 0f, 0f), 1.0f);

		Charter.launchChart("predicted", "actual", true, t1);//, t2);
	}

	public static ArrayList<XYPoint> getData(File f, int tr, int tk, int th) {
		final AminoAcidConstants aas=new AminoAcidConstants();
		final ArrayList<XYPoint> rts=new ArrayList<XYPoint>();

		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String seq=row.get("sequence");
				float actual=Float.parseFloat(row.get("actual"));
				
				float Q = getQ(seq);
				double mass=aas.getMass(seq);
				
				float predicted=(float)getCE(Q, mass);

				//if (rs==tr&&ks==tk&&hs==th) {
					rts.add(new XYPoint(predicted, actual));
				//}
			}
			
			@Override
			public void cleanup() {
			}
		};

		BlockingQueue<Map<String, String>> blockingQueue=new LinkedBlockingQueue<Map<String, String>>();
		TableParserProducer producer=new TableParserProducer(blockingQueue, f, "\t", 1);
		TableParserConsumer consumer=new TableParserConsumer(blockingQueue, muscle);

		Thread producerThread=new Thread(producer);
		Thread consumerThread=new Thread(consumer);
		producerThread.start();
		consumerThread.start();

		try {
			producerThread.join();
			consumerThread.join();
		} catch (InterruptedException ie) {
			Logger.errorLine("Percolator reading interrupted!");
			Logger.errorException(ie);
		}

		return rts;
	}
}
