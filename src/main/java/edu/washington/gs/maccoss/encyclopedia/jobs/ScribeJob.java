package edu.washington.gs.maccoss.encyclopedia.gui.framework.scribe;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import edu.washington.gs.maccoss.encyclopedia.Scribe;
import edu.washington.gs.maccoss.encyclopedia.algorithms.scribe.ScribeJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class ScribeJob extends SearchJob {
	public ScribeJob(ScribeJobData pecanData) {
		super(pecanData);
	}
	
	@Override
	public void runJob(ProgressIndicator progress) throws Exception {
		Scribe.runSearch(progress, getScribeData());
	}
	
	public ScribeJobData getScribeData() {
		return (ScribeJobData)getSearchData();
	}
	
	@Override
	public void writeToXML(Document doc, Element parentElement) {
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);
		getScribeData().writeToXML(doc, rootElement);
	}

	public static ScribeJob readFromXML(Document doc, Element rootElement) {
		if (!rootElement.getTagName().equals(ScribeJob.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+ScribeJob.class.getSimpleName()+"]");
		}

		NodeList nodes=rootElement.getChildNodes();
		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
            if (node.getNodeType() == Node.ELEMENT_NODE) {
                Element element = (Element) node;
                if (element.getTagName().equals(ScribeJobData.class.getSimpleName())) {
                	return new ScribeJob(ScribeJobData.readFromXML(doc, element));
                }
            }
		}
		throw new EncyclopediaException("Missing job data for "+rootElement.getTagName());
	}
}
