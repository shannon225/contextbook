package edu.washington.gs.maccoss.encyclopedia.jobs;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import edu.washington.gs.maccoss.encyclopedia.EncyclopediaTwo;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaTwoJobData;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class EncyclopediaTwoJob extends SearchJob {
	public EncyclopediaTwoJob(EncyclopediaTwoJobData libraryData) {
		super(libraryData);
	}
	
	@Override
	public void runJob(ProgressIndicator progress) throws Exception {
		EncyclopediaTwo.runSearch(progress, getLibraryData());
	}
	
	public EncyclopediaTwoJobData getLibraryData() {
		return (EncyclopediaTwoJobData)getSearchData();
	}
	
	@Override
	public void writeToXML(Document doc, Element parentElement) {
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);
		getLibraryData().writeToXML(doc, rootElement);
	}

	public static EncyclopediaTwoJob readFromXML(Document doc, Element rootElement) {
		if (!rootElement.getTagName().equals(EncyclopediaTwoJob.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+EncyclopediaTwoJob.class.getSimpleName()+"]");
		}

		NodeList nodes=rootElement.getChildNodes();
		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
            if (node.getNodeType() == Node.ELEMENT_NODE) {
                Element element = (Element) node;
                if (element.getTagName().equals(EncyclopediaTwoJobData.class.getSimpleName())) {
                	return new EncyclopediaTwoJob(EncyclopediaTwoJobData.readFromXML(doc, element));
                }
            }
		}
		throw new EncyclopediaException("Missing job data for "+rootElement.getTagName());
	}
}
