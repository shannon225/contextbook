package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;
import java.util.HashMap;
import java.util.Optional;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLObject;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLUtils;

public class PercolatorExecutionData implements XMLObject {
	private final File inputTSV;
	private final File fastaFile;
	private final File peptideOutputFile;
	private final File peptideDecoyFile;
	private final File proteinOutputFile;
	private final File proteinDecoyFile;
	private final SearchParameters parameters;
	private final boolean useMinMax;
	private String percolatorExecutableVersion;

	public PercolatorExecutionData(File inputTSV, File fastaFile, File peptideOutputFile, File peptideDecoyFile, File proteinOutputFile, File proteinDecoyFile, SearchParameters parameters) {
		this.inputTSV=inputTSV;
		this.fastaFile=fastaFile;
		this.peptideOutputFile=peptideOutputFile;
		this.peptideDecoyFile=peptideDecoyFile;
		this.proteinOutputFile=proteinOutputFile;
		this.proteinDecoyFile=proteinDecoyFile;
		this.parameters=parameters;
		this.useMinMax=true;
	}
	
	
	
	public PercolatorExecutionData(File inputTSV, File fastaFile, File peptideOutputFile, File peptideDecoyFile, File proteinOutputFile, File proteinDecoyFile, SearchParameters parameters, boolean useMinMax) {
		this.inputTSV=inputTSV;
		this.fastaFile=fastaFile;
		this.peptideOutputFile=peptideOutputFile;
		this.peptideDecoyFile=peptideDecoyFile;
		this.proteinOutputFile=proteinOutputFile;
		this.proteinDecoyFile=proteinDecoyFile;
		this.parameters=parameters;
		this.useMinMax=useMinMax;
	}

	@Override
	public void writeToXML(Document doc, Element parentElement) {
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);

		XMLUtils.writeTag(doc, rootElement, "inputTSV", getInputTSV().getAbsolutePath());
		XMLUtils.writeTag(doc, rootElement, "fastaFile", getFastaFile().getAbsolutePath());
		XMLUtils.writeTag(doc, rootElement, "peptideOutputFile", getPeptideOutputFile().getAbsolutePath());
		XMLUtils.writeTag(doc, rootElement, "peptideDecoyFile", getPeptideDecoyFile().getAbsolutePath());
		XMLUtils.writeTag(doc, rootElement, "proteinOutputFile", getProteinOutputFile().getAbsolutePath());
		XMLUtils.writeTag(doc, rootElement, "proteinDecoyFile", getProteinDecoyFile().getAbsolutePath());
		XMLUtils.writeTag(doc, rootElement, "useMinMax", Boolean.toString(isUseMinMax()));
		XMLUtils.writeTag(doc, rootElement, "percolatorExecutableVersion", percolatorExecutableVersion);
	}
	
	public static PercolatorExecutionData readFromXML(Document doc, Element rootElement, SearchParameters parameters) {
		if (!rootElement.getTagName().equals(PercolatorExecutionData.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+PercolatorExecutionData.class.getSimpleName()+"]");
		}
		File inputTSV=null;
		File fastaFile=null;
		File peptideOutputFile=null;
		File peptideDecoyFile=null;
		File proteinOutputFile=null;
		File proteinDecoyFile=null;
		boolean useMinMax=false;
		String percolatorExecutableVersion=null;
		
		NodeList nodes=rootElement.getChildNodes();
		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
            if (node.getNodeType() == Node.ELEMENT_NODE) {
                Element element = (Element) node;
                if ("inputTSV".equals(element.getTagName())) {
                	inputTSV=new File(element.getTextContent());
                } else if ("fastaFile".equals(element.getTagName())) {
                	fastaFile=new File(element.getTextContent());
                } else if ("peptideOutputFile".equals(element.getTagName())) {
                	peptideOutputFile=new File(element.getTextContent());
                } else if ("peptideDecoyFile".equals(element.getTagName())) {
                	peptideDecoyFile=new File(element.getTextContent());
                } else if ("proteinOutputFile".equals(element.getTagName())) {
                	proteinOutputFile=new File(element.getTextContent());
                } else if ("proteinDecoyFile".equals(element.getTagName())) {
                	proteinDecoyFile=new File(element.getTextContent());
                } else if ("useMinMax".equals(element.getTagName())) {
                	useMinMax=Boolean.parseBoolean(element.getTextContent());
                } else if ("percolatorExecutableVersion".equals(element.getTagName())) {
                	percolatorExecutableVersion=element.getTextContent();
                }
            }
		}
		
		if (inputTSV==null) throw new EncyclopediaException("Found null inputTSV in "+rootElement.getTagName());
		if (fastaFile==null) throw new EncyclopediaException("Found null fastaFile in "+rootElement.getTagName());
		if (peptideOutputFile==null) throw new EncyclopediaException("Found null peptideOutputFile in "+rootElement.getTagName());
		if (peptideDecoyFile==null) throw new EncyclopediaException("Found null peptideDecoyFile in "+rootElement.getTagName());
		if (proteinOutputFile==null) throw new EncyclopediaException("Found null proteinOutputFile in "+rootElement.getTagName());
		if (proteinDecoyFile==null) throw new EncyclopediaException("Found null proteinDecoyFile in "+rootElement.getTagName());
		if (percolatorExecutableVersion==null) throw new EncyclopediaException("Found null percolatorExecutableVersion in "+rootElement.getTagName());
		
		return new PercolatorExecutionData(inputTSV, fastaFile, peptideOutputFile, peptideDecoyFile, proteinOutputFile, proteinDecoyFile, parameters, useMinMax);
	}

	public PercolatorExecutionData getDDAVersion() {
		return new PercolatorExecutionData(inputTSV, fastaFile, peptideOutputFile, peptideDecoyFile, proteinOutputFile, proteinDecoyFile, parameters, false);
	}
	
	public boolean hasDataAvailable() {
		if (!peptideOutputFile.exists()||!peptideOutputFile.canRead()) return false;
		if (!peptideDecoyFile.exists()||!peptideDecoyFile.canRead()) return false;
		return true;
	}
	
	public File getFastaFile() {
		return fastaFile;
	}

	public File getInputTSV() {
		return inputTSV;
	}

	public File getWeightsFile(int round) {
		return new File(getPeptideOutputFile().getAbsolutePath()+"."+round+".weights");
	}

	public File getModelFile() {
		return new File(getPeptideOutputFile().getAbsolutePath()+".model");
	}

	/**
	 * This method is {@code protected} to allow access for testing.
	 *
	 * @param percolatorExecutableVersion Canonical version of Percolator parsed directly when running the actual executable
	 */
	protected void setPercolatorExecutableVersion(String percolatorExecutableVersion) {
		this.percolatorExecutableVersion = percolatorExecutableVersion;
	}

	/**
	 * @return Canonical version of Percolator parsed directly when running the actual executable
	 */
	public Optional<String> getPercolatorExecutableVersion() {
		return Optional.ofNullable(this.percolatorExecutableVersion);
	}

	public File getPeptideOutputFile() {
		return peptideOutputFile;
	}

	public File getPeptideDecoyFile() {
		return peptideDecoyFile;
	}

	public File getProteinOutputFile() {
		return proteinOutputFile;
	}

	public File getProteinDecoyFile() {
		return proteinDecoyFile;
	}
	public SearchParameters getParameters() {
		return parameters;
	}
	
	public boolean isUseMinMax() {
		return useMinMax;
	}
	
	public Optional<File> getPercolatorModelFile() {
		return parameters.getPercolatorModelFile();
	}
}
