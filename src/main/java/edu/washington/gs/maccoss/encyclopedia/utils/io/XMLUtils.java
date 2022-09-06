package edu.washington.gs.maccoss.encyclopedia.utils.io;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

public class XMLUtils {
	public static Element writeTag(Document doc, Element rootElement, String tagName, String textContent) {
		Element param=doc.createElement(tagName);
		param.setTextContent(textContent);
		rootElement.appendChild(param);
		return param;
	}
}
