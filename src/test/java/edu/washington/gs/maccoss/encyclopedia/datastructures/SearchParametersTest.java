package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.HashMap;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import junit.framework.TestCase;

public class SearchParametersTest extends TestCase {
	public void testReadParameters() throws Exception {
		HashMap<String, String> map=SearchParameterParser.getDefaultParameters();
		map.put("-runPhosphoLocalization", "true");
		map.put("-deconvoluteOverlappingWindows", "true");
		map.put("-poffset", "1000"); // definitely not default!
		
		SearchParameters params=SearchParameterParser.parseParameters(map);
		assertEquals(1000.0, params.getPrecursorOffsetPPM());
		params.savePreferences(null, null);
		
		SearchParameters readParams=SearchParameterParser.parseParameters(SearchParameters.readPreferences());
		assertEquals(1000.0, readParams.getPrecursorOffsetPPM());
		
		map.put("-poffset", "-1000"); // definitely not default!
		
		params=SearchParameterParser.parseParameters(map);
		assertEquals(-1000.0, params.getPrecursorOffsetPPM());
		params.savePreferences(null, null);
		
		readParams=SearchParameterParser.parseParameters(SearchParameters.readPreferences());
		assertEquals(-1000.0, readParams.getPrecursorOffsetPPM());
	}

	public void testEntrapmentFDR() {
		float FDR=0.01f;
		for (float i = 0.0f; i <= 2.05f; i+=0.1f) {
			System.out.println(i+"\t"+SearchParameters.getEffectivePercolatorThreshold(0.01f, i)+"\t"+(FDR * (1+i)*(1-(((1+i)-1)*FDR)))+"\t"+(FDR*(1+i)));
		}
	}
	
	public void testNonDefaultParameters() {
		SearchParameters params=SearchParameterParser.getDefaultParametersObject();
		
		HashMap<String, String> nonDefaultParameters = params.getNonDefaultParameters();
		assertTrue(nonDefaultParameters.size()==0);
		
		HashMap<String, String> map=params.toParameterMap();
		map.put("-frag", "HCD");
		map.put("-ftol", "0.5");
		map.put("-ftolunits", "AMU");
		map.put("-enzyme", "Chymotrypsin");
		map.put("-percolatorThreshold", "0.05");

		nonDefaultParameters = SearchParameterParser.parseParameters(map).getNonDefaultParameters();
//		nonDefaultParameters.forEach(new BiConsumer<String, String>() {
//			@Override
//			public void accept(String t, String u) {
//				System.out.println(t+" --> "+u);
//			}
//		});
		assertTrue(nonDefaultParameters.size()==5);
	}
	
	public void testWriteXML() throws ParserConfigurationException, TransformerException {
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
		Element rootElement = doc.createElement("root");
		doc.appendChild(rootElement);

		// add in some weird params
		PecanSearchParameters params = PecanParameterParser.getDefaultParametersObject();
		HashMap<String, String> map = params.toParameterMap();
		map.put("-frag", "HCD");
		map.put("-ftol", "0.5");
		map.put("-ftolunits", "AMU");
		map.put("-enzyme", "Chymotrypsin");
		map.put("-percolatorThreshold", "0.05");
		params = PecanParameterParser.parseParameters(map);
		
		params.writeToXML(doc, rootElement);

		TransformerFactory transformerFactory = TransformerFactory.newInstance();
		Transformer transformer = transformerFactory.newTransformer();
		transformer.setOutputProperty(OutputKeys.INDENT, "yes");

		DOMSource source = new DOMSource(doc);

		transformer.transform(source, new StreamResult(writer));
		transformer.transform(source, new StreamResult(System.out));
		
		writer.close();

		try {
	        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
	        DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
	        doc = dBuilder.parse(f);
	        doc.getDocumentElement().normalize();

			NodeList nodes=doc.getDocumentElement().getChildNodes();
			for (int i = 0; i < nodes.getLength(); i++) {
				Node node = nodes.item(i);
	            if (node.getNodeType() == Node.ELEMENT_NODE) {
	                Element element = (Element) node;
	                if (element.getTagName().equals(PecanSearchParameters.class.getSimpleName())) {
	                	PecanSearchParameters readParams=PecanSearchParameters.readFromXML(doc, element);
	        			
	        			// assert that no new or dropped params
	        			assertEquals(readParams.getSpecificParameters(params.toParameterMap()).size(), 0);
	        			assertEquals(params.getSpecificParameters(readParams.toParameterMap()).size(), 0);
	                }
	            }
			}
		} catch (IOException|SAXException e) {
			throw new EncyclopediaException("Error reading parameters!", e);
		}
	}
}
