package edu.washington.gs.maccoss.encyclopedia.algorithms.scribe;

import java.io.File;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DDASearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibrarySearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLObject;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLUtils;

public class ScribeJobData extends DDASearchJobData implements LibrarySearchJobData, XMLObject {
	public static final String LOG_FILE_SUFFIX=".log";
	public static final String OUTPUT_FILE_SUFFIX=".scribe.txt";
	public static final String DECOY_FILE_SUFFIX=".scribe.decoy.txt";
	public static final String OUTPUT_PROTEIN_FILE_SUFFIX=".scribe.protein.txt";
	public static final String DECOY_PROTEIN_FILE_SUFFIX=".scribe.protein_decoy.txt";
	public static final String FEATURE_FILE_SUFFIX=".scribe.features.txt";

	private final File fastaFile;
	private final LibraryInterface library;
	private final ScribeScoringFactory taskFactory;

	// gui
	public ScribeJobData(File diaFile, File fastaFile, LibraryInterface library, ScribeScoringFactory taskFactory) {
		this(diaFile, fastaFile, library, getPercolatorExecutionData(diaFile, fastaFile, taskFactory.getParameters()), taskFactory);
	}

	// command line
	public ScribeJobData(File diaFile, File fastaFile, LibraryInterface library, File outputFile, ScribeScoringFactory taskFactory) {
		this(diaFile, fastaFile, library, getPercolatorExecutionData(outputFile, fastaFile, taskFactory.getParameters()), taskFactory);
	}

	// internal only
	private ScribeJobData(File diaFile, File fastaFile, LibraryInterface library, PercolatorExecutionData percolatorFiles, ScribeScoringFactory taskFactory) {
		this(diaFile, null, fastaFile, library, percolatorFiles, taskFactory);
	}

	// used by testing
	public ScribeJobData(File diaFile, StripeFileInterface diaFileReader, File fastaFile, LibraryInterface library, PercolatorExecutionData percolatorFiles,
			ScribeScoringFactory taskFactory) {
		super(diaFile, diaFileReader, percolatorFiles, taskFactory.getParameters(), ProgramType.getGlobalVersion().toString());

		this.library=library;
		this.fastaFile=fastaFile;
		this.taskFactory=taskFactory;
	}
	
	@Override
	public void writeToXML(Document doc, Element parentElement) {
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);

		XMLUtils.writeTag(doc, rootElement, "fastaFile", getFastaFile().getAbsolutePath());
		XMLUtils.writeTag(doc, rootElement, "diaFile", getDiaFile().getAbsolutePath());
		XMLUtils.writeTag(doc, rootElement, "version", getVersion());
		if (library instanceof LibraryFile) {
			XMLUtils.writeTag(doc, rootElement, "library", ((LibraryFile) library).getFile().getAbsolutePath());
		}
		
		getPercolatorFiles().writeToXML(doc, rootElement);
		getParameters().writeToXML(doc, rootElement);
	}
	
	@Override
	public PercolatorExecutionData getPercolatorFiles() {
		return super.getPercolatorFiles().getDDAVersion();
	}
	
	public static ScribeJobData readFromXML(Document doc, Element rootElement) {
		if (!rootElement.getTagName().equals(ScribeJobData.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+ScribeJobData.class.getSimpleName()+"]");
		}
		File diaFile=null;
		File fastaFile=null;
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
                if (element.getTagName().equals(ScribeSearchParameters.class.getSimpleName())) {
                	readParams=ScribeSearchParameters.readFromXML(doc, element);
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
                } else if ("version".equals(element.getTagName())) {
                	version=element.getTextContent();
                } else if (element.getTagName().equals(PercolatorExecutionData.class.getSimpleName())) {
                	percolatorData=PercolatorExecutionData.readFromXML(doc, element, readParams);
                }
            }
		}

		if (fastaFile==null) throw new EncyclopediaException("Found null fastaFile in "+rootElement.getTagName());
		if (diaFile==null) throw new EncyclopediaException("Found null diaFile in "+rootElement.getTagName());
		if (library==null) throw new EncyclopediaException("Found null library in "+rootElement.getTagName());
		if (version==null) throw new EncyclopediaException("Found null version in "+rootElement.getTagName());
		if (percolatorData==null) throw new EncyclopediaException("Found null percolatorData in "+rootElement.getTagName());
		
		LibraryInterface libraryObject=BlibToLibraryConverter.getFile(library, percolatorData.getFastaFile(), readParams);

		ScribeScoringFactory factory=new ScribeScoringFactory(readParams);
		return new ScribeJobData( diaFile,  fastaFile,  libraryObject,  percolatorData,  factory);
	}
	
	@Override
	public SearchJobData updateQuantFile(File f) {
		return new ScribeJobData(f, getFastaFile(), getLibrary(), getPercolatorFiles(), getTaskFactory());
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

	public File getFastaFile() {
		return fastaFile;
	}

	public ScribeScoringFactory getTaskFactory() {
		return taskFactory;
	}

	public File getResultLibrary() {
		String absolutePath = getPrefixFromOutput(getPercolatorFiles().getPeptideOutputFile());
		return new File(absolutePath + LibraryFile.DLIB);
	}

	@Override
	public String getSearchType() {
		return "Scribe";
	}
	public LibraryInterface getLibrary() {
		return library;
	}
	
	@Override
	public String getPrimaryScoreName() {
		return taskFactory.getPrimaryScoreName();
	}
}