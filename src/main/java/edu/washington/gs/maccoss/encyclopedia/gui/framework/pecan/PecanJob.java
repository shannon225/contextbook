package edu.washington.gs.maccoss.encyclopedia.gui.framework.pecan;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import edu.washington.gs.maccoss.encyclopedia.Pecanpie;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class PecanJob extends SearchJob {
	public PecanJob(PecanJobData pecanData) {
		super(pecanData);
	}
	
	@Override
	public void runJob(ProgressIndicator progress) throws Exception {
		Pecanpie.runPie(progress, getPecanData());
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

	public static PecanJob readFromXML(Document doc, Element rootElement) {
		if (!rootElement.getTagName().equals(PecanJob.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+PecanJob.class.getSimpleName()+"]");
		}

		NodeList nodes=rootElement.getChildNodes();
		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
            if (node.getNodeType() == Node.ELEMENT_NODE) {
                Element element = (Element) node;
                if (element.getTagName().equals(PecanJobData.class.getSimpleName())) {
                	return new PecanJob(PecanJobData.readFromXML(doc, element));
                }
            }
		}
		throw new EncyclopediaException("Missing job data for "+rootElement.getTagName());
	}
}
