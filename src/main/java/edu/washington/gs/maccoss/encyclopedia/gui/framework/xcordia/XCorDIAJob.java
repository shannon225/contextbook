package edu.washington.gs.maccoss.encyclopedia.gui.framework.xcordia;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import edu.washington.gs.maccoss.encyclopedia.XCorDIA;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class XCorDIAJob extends SearchJob {
	public XCorDIAJob(JobProcessor processor, XCorDIAJobData pecanData) {
		super(processor, pecanData);
	}
	
	@Override
	public void runJob() throws Exception {
		XCorDIA.runPie(getProgressIndicator(), getXCorDIAData());
	}
	
	public XCorDIAJobData getXCorDIAData() {
		return (XCorDIAJobData)getSearchData();
	}
	
	@Override
	public void writeToXML(Document doc, Element parentElement) {
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);
		getXCorDIAData().writeToXML(doc, rootElement);
	}

	public static XCorDIAJob readFromXML(Document doc, Element rootElement, JobProcessor processor) {
		if (!rootElement.getTagName().equals(XCorDIAJob.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+XCorDIAJob.class.getSimpleName()+"]");
		}

		NodeList nodes=rootElement.getChildNodes();
		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
            if (node.getNodeType() == Node.ELEMENT_NODE) {
                Element element = (Element) node;
                if (element.getTagName().equals(XCorDIAJobData.class.getSimpleName())) {
                	return new XCorDIAJob(processor, XCorDIAJobData.readFromXML(doc, element));
                }
            }
		}
		throw new EncyclopediaException("Missing job data for "+rootElement.getTagName());
	}
}
