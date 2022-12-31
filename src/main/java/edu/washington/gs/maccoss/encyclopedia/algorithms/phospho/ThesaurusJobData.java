package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLUtils;

public class ThesaurusJobData extends EncyclopediaJobData {
	public static final String LOG_FILE_SUFFIX=".log";
	public static final String OUTPUT_FILE_SUFFIX=".thesaurus.txt";
	public static final String FEATURE_FILE_SUFFIX=".thesaurus_features.txt";
	public static final String THESAURUS_REPORT_FILE_SUFFIX=".thesaurus" + LibraryFile.ELIB;

	public ThesaurusJobData(File diaFile, LibraryInterface library, File outputFile, File fastaFile, LibraryScoringFactory taskFactory) {
		super(diaFile, getPercolatorExecutionData(outputFile, fastaFile, taskFactory.getParameters()), taskFactory.getParameters(), ProgramType.getGlobalVersion().toString(), library, taskFactory);
	}
	
	private ThesaurusJobData(File diaFile, PercolatorExecutionData percolatorFiles, SearchParameters parameters, String version, LibraryInterface library, LibraryScoringFactory taskFactory) {
		super(diaFile, percolatorFiles, parameters, version, library, taskFactory);
	}

	public static PercolatorExecutionData getPercolatorExecutionData(File referenceFileLocation, File fasta, SearchParameters parameters) {
		String prefix = getPrefix(parameters);
		return new PercolatorExecutionData(new File(getPrefixFromOutput(referenceFileLocation, parameters) + prefix+FEATURE_FILE_SUFFIX), fasta,
				new File(getPrefixFromOutput(referenceFileLocation, parameters) + prefix+OUTPUT_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation, parameters) + DECOY_FILE_SUFFIX), 
				new File(getPrefixFromOutput(referenceFileLocation, parameters) + OUTPUT_PROTEIN_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation, parameters) + DECOY_PROTEIN_FILE_SUFFIX), parameters);
	}
	
	@Override
	public void writeToXML(Document doc, Element parentElement) {
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);

		XMLUtils.writeTag(doc, rootElement, "diaFile", getDiaFile().getAbsolutePath());
		XMLUtils.writeTag(doc, rootElement, "version", getVersion());
		if (getLibrary() instanceof LibraryFile) {
			XMLUtils.writeTag(doc, rootElement, "library", ((LibraryFile) getLibrary()).getFile().getAbsolutePath());
		}

		getPercolatorFiles().writeToXML(doc, rootElement);
		getParameters().writeToXML(doc, rootElement);
	}


	public static ThesaurusJobData readFromXML(Document doc, Element rootElement) {
		if (!rootElement.getTagName().equals(ThesaurusJobData.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+ThesaurusJobData.class.getSimpleName()+"]");
		}
		File diaFile=null;
		File library=null;
		String version=null;
		PercolatorExecutionData percolatorData=null;
		ThesaurusSearchParameters readParams=null;

		NodeList nodes=rootElement.getChildNodes();

		// read params first
		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
            if (node.getNodeType() == Node.ELEMENT_NODE) {
                Element element = (Element) node;
                if (element.getTagName().equals(SearchParameters.class.getSimpleName())) {
                	readParams=ThesaurusSearchParameters.readFromXML(doc, element);
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
		return new ThesaurusJobData(diaFile,  percolatorData, readParams,  version, libraryObject, factory);
	}

	@Override
	public SearchJobData updateQuantFile(File f) {
		return new ThesaurusJobData(f, getPercolatorFiles(), getParameters(), getVersion(), getLibrary(), getTaskFactory());
	}

	public ThesaurusJobData updateTaskFactory(LibraryScoringFactory taskFactory) {
		return new ThesaurusJobData(getDiaFile(), getPercolatorFiles(), getParameters(), getVersion(), getLibrary(), taskFactory);
	}
	
	public File getLocalizationFile() {
		String prefix = getPrefix(getParameters());
		String absolutePath = getPrefixFromOutput(getPercolatorFiles().getPeptideOutputFile(), getParameters());
		return new File(absolutePath+prefix+".localizations.txt");
	}

	public File getResultLibrary() {
		String prefix = getPrefix(getParameters());
		String absolutePath = getPrefixFromOutput(getPercolatorFiles().getPeptideOutputFile(), getParameters());
		return new File(absolutePath+prefix+THESAURUS_REPORT_FILE_SUFFIX);
	}

	static String getPrefixFromOutput(File outputFile, SearchParameters parameters) {
		final String absolutePath = outputFile.getAbsolutePath();

		String prefix = getPrefix(parameters);
		if (absolutePath.endsWith(prefix+OUTPUT_FILE_SUFFIX)) {
			return absolutePath.substring(0, absolutePath.length() - (prefix.length()+OUTPUT_FILE_SUFFIX.length()));
		} else {
			return absolutePath;
		}
	}
	
	@Override
	public String getSearchType() {
		return "Thesaurus";
	}

	private static String getPrefix(SearchParameters parameters) {
		if (!parameters.getLocalizingModification().isPresent()) {
			throw new EncyclopediaException("Sorry, Thesaurus requires you to specify a localizing PTM!");
		}
		String prefix="."+parameters.getLocalizingModification().get().getShortname();
		return prefix;
	}
}