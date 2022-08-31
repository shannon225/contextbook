package edu.washington.gs.maccoss.encyclopedia.gui.framework.library;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import edu.washington.gs.maccoss.encyclopedia.Thesaurus;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ThesaurusJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class ThesaurusJob extends SearchJob {
	public ThesaurusJob(ThesaurusJobData libraryData) {
		super(libraryData);
	}
	
	@Override
	public void runJob(ProgressIndicator progress) throws Exception {
		Thesaurus.runSearch(progress, getLibraryData());
	}
	
	public ThesaurusJobData getLibraryData() {
		return (ThesaurusJobData)getSearchData();
	}
	
	@Override
	public void writeToXML(Document doc, Element parentElement) {
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);
		getLibraryData().writeToXML(doc, rootElement);
	}

	public static ThesaurusJob readFromXML(Document doc, Element rootElement) {
		if (!rootElement.getTagName().equals(ThesaurusJob.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+ThesaurusJob.class.getSimpleName()+"]");
		}

		NodeList nodes=rootElement.getChildNodes();
		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
            if (node.getNodeType() == Node.ELEMENT_NODE) {
                Element element = (Element) node;
                if (element.getTagName().equals(ThesaurusJobData.class.getSimpleName())) {
                	return new ThesaurusJob(ThesaurusJobData.readFromXML(doc, element));
                }
            }
		}
		throw new EncyclopediaException("Missing job data for "+rootElement.getTagName());
	}
}
