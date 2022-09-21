package edu.washington.gs.maccoss.encyclopedia.utils.io;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

public interface XMLObject {
	public void writeToXML(Document doc, Element parentElement);
}
