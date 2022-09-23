package edu.washington.gs.maccoss.encyclopedia.jobs;

import edu.washington.gs.maccoss.encyclopedia.SearchToBLIB;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLObject;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import java.io.File;
import java.util.ArrayList;

public class AlignmentOnlyLibraryJob implements WorkerJob, XMLObject {
	private final File destFile;
	private final JobProcessor processor;

	public AlignmentOnlyLibraryJob(File destFile, JobProcessor processor) {
		this.destFile = destFile;
		this.processor = processor;
	}

	@Override
	public void runJob(ProgressIndicator progress) throws Exception {
		ArrayList<SearchJobData> jobData=new ArrayList<SearchJobData>();
		for (WorkerJob job : processor.getQueue()) {
			if (job instanceof SearchJob) {
				jobData.add(((SearchJob)job).getSearchData());
			}
		}

		Logger.logLine("Found "+jobData.size()+" jobs in the queue to combine...");

		if (!jobData.isEmpty()) {
			final SearchJobData reprJob = jobData.iterator().next();
			Logger.logLine("Extracting representative search parameters from job " + reprJob.getDiaFileReader().getOriginalFileName());

			SearchToBLIB.convert(
					progress,
					jobData,
					destFile,
					SearchToBLIB.OutputFormat.ALIB,
					true, // Must be true for ALIB export
					reprJob.getParameters()
			);
		}
	}

	@Override
	public String getJobTitle() {
		return "Export alignment-only library " + destFile.getAbsolutePath();
	}

	@Override
	public void writeToXML(Document doc, Element parentElement) {
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);

		XMLUtils.writeTag(doc, rootElement, "elibFile", destFile.getAbsolutePath());
	}

	public static AlignmentOnlyLibraryJob readFromXML(Document doc, Element rootElement, JobProcessor processor) {
		if (!rootElement.getTagName().equals(AlignmentOnlyLibraryJob.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+AlignmentOnlyLibraryJob.class.getSimpleName()+"]");
		}
		File elibFile=null;

		NodeList nodes=rootElement.getChildNodes();

		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
			if (node.getNodeType() == Node.ELEMENT_NODE) {
				Element element = (Element) node;
				if ("elibFile".equals(element.getTagName())) {
					elibFile=new File(element.getTextContent());
				}
			}
		}

		if (elibFile==null) throw new EncyclopediaException("Found null elibFile in "+rootElement.getTagName());
		return new AlignmentOnlyLibraryJob(elibFile, processor);
	}
}
