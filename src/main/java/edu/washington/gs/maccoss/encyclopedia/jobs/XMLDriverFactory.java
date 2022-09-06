package edu.washington.gs.maccoss.encyclopedia.jobs;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLObject;

public class XMLDriverFactory {
	public static String DRIVER_XML_EXTENSION = ".encxml";
	
	public static void writeXML(ArrayList<WorkerJob> jobs, File file) {
		PrintWriter writer=null;
		try {
			writer=new PrintWriter(file, "UTF-8");

			DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
			DocumentBuilder docBuilder = docFactory.newDocumentBuilder();

			// root elements
			Document doc = docBuilder.newDocument();
			Element rootElement = doc.createElement("EncyclopeDIA");
			doc.appendChild(rootElement);

			for (WorkerJob job : jobs) {
				if (job instanceof XMLObject) {
					((XMLObject)job).writeToXML(doc, rootElement);
				}
			}

			TransformerFactory transformerFactory = TransformerFactory.newInstance();
			Transformer transformer = transformerFactory.newTransformer();
			transformer.setOutputProperty(OutputKeys.INDENT, "yes");

			DOMSource source = new DOMSource(doc);

			transformer.transform(source, new StreamResult(writer));
			transformer.transform(source, new StreamResult(System.out));
			
			
		} catch (Exception e) {
			throw new EncyclopediaException("Error writing parameters!", e);
		} finally {
			if (writer!=null) {
				writer.close();
			}
		}
	}
	
	public static ArrayList<WorkerJob> readXML(File file) {
		try {
			ArrayList<WorkerJob> jobs=new ArrayList<WorkerJob>();
			
			DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
			DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
			Document doc = dBuilder.parse(file);
			doc.getDocumentElement().normalize();

			NodeList nodes = doc.getDocumentElement().getChildNodes();
			for (int i = 0; i < nodes.getLength(); i++) {
				Node node = nodes.item(i);
				if (node.getNodeType() == Node.ELEMENT_NODE) {
					Element element = (Element) node;
					
					if (EncyclopediaJob.class.getSimpleName().equals(element.getTagName())) {
						jobs.add(EncyclopediaJob.readFromXML(doc, element));
					} else if (ThesaurusJob.class.getSimpleName().equals(element.getTagName())) {
						jobs.add(ThesaurusJob.readFromXML(doc, element));
					} else if (PecanJob.class.getSimpleName().equals(element.getTagName())) {
						jobs.add(PecanJob.readFromXML(doc, element));
					} else if (XCorDIAJob.class.getSimpleName().equals(element.getTagName())) {
						jobs.add(XCorDIAJob.readFromXML(doc, element));
					} else if (ScribeJob.class.getSimpleName().equals(element.getTagName())) {
						jobs.add(ScribeJob.readFromXML(doc, element));
					}
				}
			}
			
			return jobs;
		} catch (Exception e) {
			throw new EncyclopediaException("Error reading parameters!", e);
		}
	}
}
