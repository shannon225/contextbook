package edu.washington.gs.maccoss.encyclopedia.jobs;

import edu.washington.gs.maccoss.encyclopedia.SearchToBLIB;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import java.io.File;

public class AlignmentOnlyLibraryJob extends LibExportJob {
	public AlignmentOnlyLibraryJob(File destFile, JobProcessor processor) {
		super(
				destFile,
				SearchToBLIB.OutputFormat.ALIB,
				true, // alignBetweenFiles MUST be true for ALIB export
				processor
		);
	}

	@Override
	public String getJobTitle() {
		return "Export alignment-only library " + destFile.getAbsolutePath();
	}

	@Override
	public Element writeToXML(Document doc, Element parentElement) {
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);

		XMLUtils.writeTag(doc, rootElement, "elibFile", destFile.getAbsolutePath());
		return rootElement;
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
