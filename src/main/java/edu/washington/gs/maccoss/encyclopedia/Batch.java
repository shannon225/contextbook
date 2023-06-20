package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.parameters.InstrumentSpecificSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.jobs.WorkerJob;
import edu.washington.gs.maccoss.encyclopedia.jobs.XMLDriverFactory;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class Batch {

	public static void main(String[] args) {
		HashMap<String, String> arguments=CommandLineParser.parseArguments(args);
		arguments=InstrumentSpecificSearchParameters.checkParameters(arguments);

		if (arguments.containsKey("-h")||arguments.containsKey("-help")||arguments.containsKey("--help")) {
			Logger.logLine("EncyclopeDIA Batch Help");
			Logger.timelessLogLine("EncyclopeDIA Batch allows to driving EncyclopeDIA with an XML file (*.encxml).");
			Logger.timelessLogLine("These files can be created by the GUI by saving an XML driver file of a queue.");
		}
		
		File f=new File(arguments.get("-batch"));
		if (f.exists()&&f.canRead()) {
			ArrayList<WorkerJob> jobs=XMLDriverFactory.readXML(f);
			ProgressIndicator progress=new EmptyProgressIndicator(false);
			
			for (WorkerJob job : jobs) {
				try {
					job.runJob(progress);
				} catch (Exception e) {
					throw new EncyclopediaException("Encountered unexpected exception running "+job.getJobTitle(), e);
				}
			}
			
		} else {
			Logger.errorLine("You are required to specify an driver XML file ("+XMLDriverFactory.DRIVER_XML_EXTENSION+")");
		}
	}

}
