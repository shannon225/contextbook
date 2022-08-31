package edu.washington.gs.maccoss.encyclopedia.gui.framework.pecan;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import edu.washington.gs.maccoss.encyclopedia.Pecanpie;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class PecanJob extends SearchJob {
	public PecanJob(JobProcessor processor, PecanJobData pecanData) {
		super(processor, pecanData);
	}
	
	@Override
	public void runJob() throws Exception {
		Pecanpie.runPie(getProgressIndicator(), getPecanData());
	}
	
	public PecanJobData getPecanData() {
		return (PecanJobData)getSearchData();
	}
	
	@Override
	public void writeToXML(Document doc, Element parentElement) {
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);
		getPecanData().writeToXML(doc, rootElement);
	}

	public static PecanJob readFromXML(Document doc, Element rootElement, JobProcessor processor) {
		if (!rootElement.getTagName().equals(PecanJob.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+PecanJob.class.getSimpleName()+"]");
		}

		NodeList nodes=rootElement.getChildNodes();
		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
            if (node.getNodeType() == Node.ELEMENT_NODE) {
                Element element = (Element) node;
                if (element.getTagName().equals(PecanJobData.class.getSimpleName())) {
                	return new PecanJob(processor, PecanJobData.readFromXML(doc, element));
                }
            }
		}
		throw new EncyclopediaException("Missing job data for "+rootElement.getTagName());
	}
}
