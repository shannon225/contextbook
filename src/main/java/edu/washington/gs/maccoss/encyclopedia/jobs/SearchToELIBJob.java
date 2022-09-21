package edu.washington.gs.maccoss.encyclopedia.jobs;

import java.io.File;
import java.util.ArrayList;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import edu.washington.gs.maccoss.encyclopedia.SearchToBLIB;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLObject;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class SearchToELIBJob implements WorkerJob, XMLObject {
	private final File elibFile;
	private final boolean alignBetweenFiles;
	JobProcessor processor;

	public SearchToELIBJob(File elibFile, boolean alignBetweenFiles, JobProcessor processor) {
		this.processor=processor;
		this.elibFile=elibFile;
		this.alignBetweenFiles=alignBetweenFiles;
	}
	
	@Override
	public String getJobTitle() {
		return "Write Library "+elibFile.getName();
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
		SearchToBLIB.convert(progress, jobData, elibFile, false, alignBetweenFiles);
	}
	
	@Override
	public void writeToXML(Document doc, Element parentElement) {
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);

		XMLUtils.writeTag(doc, rootElement, "elibFile", elibFile.getAbsolutePath());
		XMLUtils.writeTag(doc, rootElement, "alignBetweenFiles", Boolean.toString(alignBetweenFiles));
	}
	
	public static SearchToELIBJob readFromXML(Document doc, Element rootElement, JobProcessor processor) {
		if (!rootElement.getTagName().equals(SearchToELIBJob.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+SearchToELIBJob.class.getSimpleName()+"]");
		}
		File elibFile=null;
		boolean alignBetweenFiles=false;
		
		NodeList nodes=rootElement.getChildNodes();

		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
            if (node.getNodeType() == Node.ELEMENT_NODE) {
                Element element = (Element) node;
                if ("elibFile".equals(element.getTagName())) {
                	elibFile=new File(element.getTextContent());
                } else if ("alignBetweenFiles".equals(element.getTagName())) {
                	alignBetweenFiles=Boolean.parseBoolean(element.getTextContent());
                }
            }
		}
		
		if (elibFile==null) throw new EncyclopediaException("Found null elibFile in "+rootElement.getTagName());
		return new SearchToELIBJob(elibFile, alignBetweenFiles, processor);
	}
}
