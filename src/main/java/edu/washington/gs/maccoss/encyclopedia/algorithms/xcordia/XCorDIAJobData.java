package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.QuantitativeSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLObject;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLUtils;

public class XCorDIAJobData extends QuantitativeSearchJobData implements XMLObject {
	public static final String OUTPUT_FILE_SUFFIX=".xcordia.txt";
	public static final String DECOY_FILE_SUFFIX=".xcordia.decoy.txt";
	public static final String OUTPUT_PROTEIN_FILE_SUFFIX=".xcordia.protein.txt";
	public static final String DECOY_PROTEIN_FILE_SUFFIX=".xcordia.protein_decoy.txt";
	public static final String FEATURE_FILE_SUFFIX=".features.txt";

	private final Optional<ArrayList<FastaPeptideEntry>> targetList;
	private final File fastaFile;
	private final Optional<LibraryInterface> library;
	private final XCorDIAOneScoringFactory taskFactory;

	// gui
	public XCorDIAJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, Optional<LibraryInterface> library, File diaFile, File fastaFile, XCorDIAOneScoringFactory taskFactory) {
		this(targetList, library, diaFile, fastaFile, getPercolatorExecutionData(diaFile, fastaFile, taskFactory.getParameters()), taskFactory);
	}

	// command line
	public XCorDIAJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, Optional<LibraryInterface> library, File diaFile, File fastaFile, File outputFile, XCorDIAOneScoringFactory taskFactory) {
		this(targetList, library, diaFile, fastaFile, getPercolatorExecutionData(outputFile, fastaFile, taskFactory.getParameters()), taskFactory);
	}

	// internal only
	protected XCorDIAJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, Optional<LibraryInterface> library, File diaFile, File fastaFile, PercolatorExecutionData percolatorFiles, XCorDIAOneScoringFactory taskFactory) {
		this(targetList, library, diaFile, null, fastaFile, percolatorFiles, taskFactory);
	}

	// used by testing
	public XCorDIAJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, Optional<LibraryInterface> library, File diaFile, StripeFileInterface diaFileReader, File fastaFile, PercolatorExecutionData percolatorFiles,
			XCorDIAOneScoringFactory taskFactory) {
		super(diaFile, diaFileReader, percolatorFiles, taskFactory.getParameters(), ProgramType.getGlobalVersion().toString());

		this.targetList=targetList;
		this.library=library;
		this.fastaFile=fastaFile;
		this.taskFactory=taskFactory;
	}
	
	@Override
	public void writeToXML(Document doc, Element parentElement) {
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);

		XMLUtils.writeTag(doc, rootElement, "diaFile", getDiaFile().getAbsolutePath());
		XMLUtils.writeTag(doc, rootElement, "fastaFile", getFastaFile().getAbsolutePath());
		if (library.isPresent()&&library.get() instanceof LibraryFile) {
			XMLUtils.writeTag(doc, rootElement, "library", ((LibraryFile) library.get()).getFile().getAbsolutePath());
		}
		XMLUtils.writeTag(doc, rootElement, "version", getVersion());
		
		getPercolatorFiles().writeToXML(doc, rootElement);
		getParameters().writeToXML(doc, rootElement);
		
		if (getTargetList().isPresent()) {
			Element targetElement=doc.createElement("targetList");
			rootElement.appendChild(targetElement);
			for (FastaPeptideEntry peptide : getTargetList().get()) {
				peptide.writeToXML(doc, parentElement);
			}
		}
	}

	
	public static XCorDIAJobData readFromXML(Document doc, Element rootElement) {
		if (!rootElement.getTagName().equals(XCorDIAJobData.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+XCorDIAJobData.class.getSimpleName()+"]");
		}
		File diaFile=null;
		File fastaFile=null;
		File library=null;
		PercolatorExecutionData percolatorData=null;
		PecanSearchParameters readParams=null;
		ArrayList<FastaPeptideEntry> targetList=new ArrayList<FastaPeptideEntry>();
		
		NodeList nodes=rootElement.getChildNodes();

		// read params first
		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
            if (node.getNodeType() == Node.ELEMENT_NODE) {
                Element element = (Element) node;
                if (element.getTagName().equals(PecanSearchParameters.class.getSimpleName())) {
                	readParams=PecanSearchParameters.readFromXML(doc, element);
                }
            }
		}
		if (readParams==null) throw new EncyclopediaException("Found null readParams in "+rootElement.getTagName());
		
		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
            if (node.getNodeType() == Node.ELEMENT_NODE) {
                Element element = (Element) node;
                if ("diaFile".equals(element.getTagName())) {
                	diaFile=new File(element.getTextContent());
                } else if ("library".equals(element.getTagName())) {
                	library=new File(element.getTextContent());
                } else if ("fastaFile".equals(element.getTagName())) {
                	fastaFile=new File(element.getTextContent());
                } else if ("targetList".equals(element.getTagName())) {
                	NodeList targets=element.getChildNodes();
            		for (int j = 0; j < targets.getLength(); j++) {
            			Node targetNode = targets.item(j);
                        if (targetNode.getNodeType() == Node.ELEMENT_NODE) {
                            Element targetNodeElement = (Element) targetNode;
                            targetList.add(FastaPeptideEntry.readFromXML(doc, targetNodeElement));
                        }
            		}
                } else if (element.getTagName().equals(PercolatorExecutionData.class.getSimpleName())) {
                	percolatorData=PercolatorExecutionData.readFromXML(doc, element, readParams);
                }
            }
		}
		
		if (diaFile==null) throw new EncyclopediaException("Found null diaFile in "+rootElement.getTagName());
		if (fastaFile==null) throw new EncyclopediaException("Found null fastaFile in "+rootElement.getTagName());
		if (percolatorData==null) throw new EncyclopediaException("Found null percolatorData in "+rootElement.getTagName());
		
		if (targetList.size()==0) targetList=null;

		LibraryInterface libraryObject=null;
		if (library!=null) {
			libraryObject=BlibToLibraryConverter.getFile(library, percolatorData.getFastaFile(), readParams);
		}
		XCorDIAOneScoringFactory factory=new XCorDIAOneScoringFactory(readParams);
		return new XCorDIAJobData(Optional.ofNullable(targetList), Optional.ofNullable(libraryObject), diaFile, fastaFile, percolatorData, factory);
	}
	
	@Override
	public SearchJobData updateQuantFile(File f) {
		return new XCorDIAJobData(getTargetList(), getLibrary(), f, getFastaFile(), getPercolatorFiles(), getTaskFactory());
	}

	public static PercolatorExecutionData getPercolatorExecutionData(File referenceFileLocation, File fastaFile, SearchParameters parameters) {
		return new PercolatorExecutionData(new File(getPrefixFromOutput(referenceFileLocation) + FEATURE_FILE_SUFFIX), fastaFile, 
				new File(getPrefixFromOutput(referenceFileLocation) + OUTPUT_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation) + DECOY_FILE_SUFFIX), 
				new File(getPrefixFromOutput(referenceFileLocation) + OUTPUT_PROTEIN_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation) + DECOY_PROTEIN_FILE_SUFFIX), parameters);
	}

	protected static String getPrefixFromOutput(File outputFile) {
		final String absolutePath = outputFile.getAbsolutePath();

		if (absolutePath.endsWith(OUTPUT_FILE_SUFFIX)) {
			return absolutePath.substring(0, absolutePath.length() - OUTPUT_FILE_SUFFIX.length());
		} else {
			return absolutePath;
		}
	}

	public Optional<ArrayList<FastaPeptideEntry>> getTargetList() {
		return targetList;
	}

	public File getFastaFile() {
		return fastaFile;
	}

	public XCorDIAOneScoringFactory getTaskFactory() {
		return taskFactory;
	}

	public File getResultLibrary() {
		String absolutePath = getPrefixFromOutput(getPercolatorFiles().getPeptideOutputFile());
		return new File(absolutePath + LibraryFile.ELIB);
	}

	@Override
	public String getSearchType() {
		return "XCorDIA";
	}
	public Optional<LibraryInterface> getLibrary() {
		return library;
	}
	
	@Override
	public String getPrimaryScoreName() {
		return taskFactory.getPrimaryScoreName();
	}
}