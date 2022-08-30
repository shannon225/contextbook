package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibrarySearchJobData;
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

public class EncyclopediaJobData extends QuantitativeSearchJobData implements LibrarySearchJobData, XMLObject {
	public static final String LOG_FILE_SUFFIX=".log";
	public static final String DECOY_PROTEIN_FILE_SUFFIX=".encyclopedia.protein_decoy.txt";
	public static final String OUTPUT_PROTEIN_FILE_SUFFIX=".encyclopedia.protein.txt";
	public static final String DECOY_FILE_SUFFIX=".encyclopedia.decoy.txt";
	public static final String OUTPUT_FILE_SUFFIX=".encyclopedia.txt";
	public static final String FEATURE_FILE_SUFFIX=".features.txt";

	private final LibraryInterface library;
	private final LibraryScoringFactory taskFactory;

	public EncyclopediaJobData(File diaFile, File fastaFile, LibraryInterface library, LibraryScoringFactory taskFactory) {
		this(diaFile, null, getPercolatorExecutionData(diaFile, fastaFile, taskFactory.getParameters()), taskFactory.getParameters(), ProgramType.getGlobalVersion().toString(), library, taskFactory);
	}

	public EncyclopediaJobData(File diaFile, File fastaFile, LibraryInterface library, File outputFile, LibraryScoringFactory taskFactory) {
		this(diaFile, null, getPercolatorExecutionData(outputFile, fastaFile, taskFactory.getParameters()), taskFactory.getParameters(), ProgramType.getGlobalVersion().toString(), library, taskFactory);
	}

	public EncyclopediaJobData(File diaFile, PercolatorExecutionData percolatorFiles, SearchParameters parameters, String version, LibraryInterface library, LibraryScoringFactory taskFactory) {
		this(diaFile, null, percolatorFiles, parameters, version, library, taskFactory);
	}

	public EncyclopediaJobData(File diaFile, StripeFileInterface diaFileReader, PercolatorExecutionData percolatorFiles, SearchParameters parameters, String version, LibraryInterface library, LibraryScoringFactory taskFactory) {
		super(diaFile, diaFileReader, percolatorFiles, parameters, version);

		this.library = library;
		this.taskFactory = taskFactory;
	}
	
	@Override
	public void writeToXML(Document doc, Element parentElement) {
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);

		XMLUtils.writeTag(doc, rootElement, "diaFile", getDiaFile().getAbsolutePath());
		XMLUtils.writeTag(doc, rootElement, "version", getVersion());
		if (library instanceof LibraryFile) {
			XMLUtils.writeTag(doc, rootElement, "library", ((LibraryFile) library).getFile().getAbsolutePath());
		}
		
		getPercolatorFiles().writeToXML(doc, rootElement);
		getParameters().writeToXML(doc, rootElement);
	}

	
	public static EncyclopediaJobData readFromXML(Document doc, Element rootElement) {
		if (!rootElement.getTagName().equals(EncyclopediaJobData.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+EncyclopediaJobData.class.getSimpleName()+"]");
		}
		File diaFile=null;
		File library=null;
		String version=null;
		PercolatorExecutionData percolatorData=null;
		SearchParameters readParams=null;
		
		NodeList nodes=rootElement.getChildNodes();

		// read params first
		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
            if (node.getNodeType() == Node.ELEMENT_NODE) {
                Element element = (Element) node;
                if (element.getTagName().equals(SearchParameters.class.getSimpleName())) {
                	readParams=SearchParameters.readFromXML(doc, element);
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
                } else if ("version".equals(element.getTagName())) {
                	version=element.getTextContent();
                } else if (element.getTagName().equals(PercolatorExecutionData.class.getSimpleName())) {
                	percolatorData=PercolatorExecutionData.readFromXML(doc, element, readParams);
                }
            }
		}
		
		if (diaFile==null) throw new EncyclopediaException("Found null diaFile in "+rootElement.getTagName());
		if (library==null) throw new EncyclopediaException("Found null library in "+rootElement.getTagName());
		if (version==null) throw new EncyclopediaException("Found null version in "+rootElement.getTagName());
		if (percolatorData==null) throw new EncyclopediaException("Found null percolatorData in "+rootElement.getTagName());
		
		LibraryInterface libraryObject=BlibToLibraryConverter.getFile(library, percolatorData.getFastaFile(), readParams);

		LibraryScoringFactory factory=EncyclopediaScoringFactory.getDefaultScoringFactory(readParams);
		return new EncyclopediaJobData(diaFile,  percolatorData, readParams,  version, libraryObject, factory);
	}
	
	@Override
	public SearchJobData updateQuantFile(File f) {
		return new EncyclopediaJobData(f, getPercolatorFiles(), getParameters(), getVersion(), getLibrary(), getTaskFactory());
	}

	protected static PercolatorExecutionData getPercolatorExecutionData(File referenceFileLocation, File fastaFile, SearchParameters parameters) {
		return new PercolatorExecutionData(new File(getPrefixFromOutput(referenceFileLocation) + FEATURE_FILE_SUFFIX), fastaFile,
				new File(getPrefixFromOutput(referenceFileLocation) + OUTPUT_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation) + DECOY_FILE_SUFFIX), 
				new File(getPrefixFromOutput(referenceFileLocation) + OUTPUT_PROTEIN_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation) + DECOY_PROTEIN_FILE_SUFFIX), parameters);
	}

	static String getPrefixFromOutput(File outputFile) {
		final String absolutePath = outputFile.getAbsolutePath();

		if (absolutePath.endsWith(OUTPUT_FILE_SUFFIX)) {
			return absolutePath.substring(0, absolutePath.length() - OUTPUT_FILE_SUFFIX.length());
		} else {
			return absolutePath;
		}
	}

	public EncyclopediaJobData updateTaskFactory(LibraryScoringFactory taskFactory) {
		return new EncyclopediaJobData(getDiaFile(), diaFileReader, getPercolatorFiles(), getParameters(), getVersion(), getLibrary(), taskFactory);
	}

	public LibraryInterface getLibrary() {
		return library;
	}

	public LibraryScoringFactory getTaskFactory() {
		return taskFactory;
	}

	public File getResultLibrary() {
		String absolutePath = getPrefixFromOutput(getPercolatorFiles().getPeptideOutputFile());
		return new File(absolutePath + LibraryFile.ELIB);
	}

	@Override
	public String getSearchType() {
		return "EncyclopeDIA";
	}
	
	@Override
	public String getPrimaryScoreName() {
		return taskFactory.getPrimaryScoreName();
	}
}