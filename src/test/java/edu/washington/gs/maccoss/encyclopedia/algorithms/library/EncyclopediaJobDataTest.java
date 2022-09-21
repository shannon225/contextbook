package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;
import java.io.PrintWriter;
import java.util.HashMap;

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

import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.algorithms.SearchTestSupport;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import junit.framework.TestCase;

public class EncyclopediaJobDataTest extends TestCase {
	public void testGetOutputAbsolutePathPrefix() {
		String absolutePath="/Users/searleb/Documents/projects/encyclopedia/encyclopedia/121115_BCS_HeLa_24mz_400_1000.dia.encyclopedia.txt";
		assertEquals("/Users/searleb/Documents/projects/encyclopedia/encyclopedia/121115_BCS_HeLa_24mz_400_1000.dia", EncyclopediaJobData.getPrefixFromOutput(new File(absolutePath)));

		absolutePath="/Users/searleb/Documents/projects/encyclopedia/encyclopedia/121115_BCS_HeLa_24mz_400_1000.dia.encyclopedia";
		assertEquals("/Users/searleb/Documents/projects/encyclopedia/encyclopedia/121115_BCS_HeLa_24mz_400_1000.dia.encyclopedia", EncyclopediaJobData.getPrefixFromOutput(new File(absolutePath)));
	}
	
	public void testXML() throws Exception {
		SearchParameters params = SearchParameterParser.getDefaultParametersObject();
		HashMap<String, String> map = params.toParameterMap();
		map.put("-frag", "HCD");
		map.put("-ftol", "0.5");
		map.put("-ftolunits", "AMU");
		map.put("-enzyme", "Chymotrypsin");
		map.put("-percolatorThreshold", "0.05");
		params = SearchParameterParser.parseParameters(map);
		
		File diaFile=SearchTestSupport.getTestFullWindowFile();
		File fastaFile=SearchTestSupport.getTestFastaFile();
		LibraryFile library=SearchTestSupport.getLibrary();
		File outputFile=new File(diaFile.getAbsolutePath()+LibraryFile.ELIB);
		
		LibraryScoringFactory factory=EncyclopediaScoringFactory.getDefaultScoringFactory(params);
		
		String version=ProgramType.getGlobalVersion().toString();
		PercolatorExecutionData percolatorData=EncyclopediaJobData.getPercolatorExecutionData(outputFile, fastaFile, factory.getParameters());
		EncyclopediaJobData job=new EncyclopediaJobData(diaFile, percolatorData, params, version, library, factory);

		File f;
		PrintWriter writer;
		try {
			f=File.createTempFile("parameters", ".xml");
			f.deleteOnExit();
			writer=new PrintWriter(f, "UTF-8");
		} catch (Exception e) {
			throw new EncyclopediaException("Error writing parameters!", e);
		}

		DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
		DocumentBuilder docBuilder = docFactory.newDocumentBuilder();

		// root elements
		Document doc = docBuilder.newDocument();
		Element rootElement = doc.createElement("EncyclopeDIA");
		doc.appendChild(rootElement);
		
		job.writeToXML(doc, rootElement);

		TransformerFactory transformerFactory = TransformerFactory.newInstance();
		Transformer transformer = transformerFactory.newTransformer();
		transformer.setOutputProperty(OutputKeys.INDENT, "yes");

		DOMSource source = new DOMSource(doc);

		transformer.transform(source, new StreamResult(writer));
		transformer.transform(source, new StreamResult(System.out));
		
		writer.close();
		
        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
        doc = dBuilder.parse(f);
        doc.getDocumentElement().normalize();

		NodeList nodes=doc.getDocumentElement().getChildNodes();
		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
            if (node.getNodeType() == Node.ELEMENT_NODE) {
                Element element = (Element) node;
                assertEquals(EncyclopediaJobData.class.getSimpleName(), element.getTagName());
            	EncyclopediaJobData data=EncyclopediaJobData.readFromXML(doc, element);
                assertNotNull(data);
            }
		}
	}
}
