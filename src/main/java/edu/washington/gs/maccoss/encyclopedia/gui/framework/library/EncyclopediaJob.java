package edu.washington.gs.maccoss.encyclopedia.gui.framework.library;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import edu.washington.gs.maccoss.encyclopedia.Encyclopedia;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class EncyclopediaJob extends SearchJob {
	public EncyclopediaJob(EncyclopediaJobData libraryData) {
		super(libraryData);
	}
	
	@Override
	public void runJob(ProgressIndicator progress) throws Exception {
		Encyclopedia.runSearch(progress, getLibraryData());
	}
	
	public EncyclopediaJobData getLibraryData() {
		return (EncyclopediaJobData)getSearchData();
	}
	
	@Override
	public void writeToXML(Document doc, Element parentElement) {
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);
		getLibraryData().writeToXML(doc, rootElement);
	}

	public static EncyclopediaJob readFromXML(Document doc, Element rootElement) {
		if (!rootElement.getTagName().equals(EncyclopediaJob.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+EncyclopediaJob.class.getSimpleName()+"]");
		}

		NodeList nodes=rootElement.getChildNodes();
		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
            if (node.getNodeType() == Node.ELEMENT_NODE) {
                Element element = (Element) node;
                if (element.getTagName().equals(EncyclopediaJobData.class.getSimpleName())) {
                	return new EncyclopediaJob(EncyclopediaJobData.readFromXML(doc, element));
                }
            }
		}
		throw new EncyclopediaException("Missing job data for "+rootElement.getTagName());
	}
}
