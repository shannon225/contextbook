package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Optional;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.zip.DataFormatException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.BackgroundFrequencyCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.BackgroundFrequencyInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.VariantXCorDIAOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLUtils;

public class VariantXCorDIAJobData extends XCorDIAJobData {
	// GUI
	public VariantXCorDIAJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, Optional<LibraryInterface> library, File diaFile, File fastaFile, XCorDIAOneScoringFactory taskFactory) {
		super(targetList, library, diaFile, fastaFile, taskFactory);
	}

	// command line
	public VariantXCorDIAJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, Optional<LibraryInterface> library, File diaFile, File fastaFile, File outputFile, XCorDIAOneScoringFactory taskFactory) {
		super(targetList, library, diaFile, fastaFile, outputFile, taskFactory);
	}
	// internal
	public VariantXCorDIAJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, Optional<LibraryInterface> library, File diaFile, File fastaFile, PercolatorExecutionData percolatorFiles,
			XCorDIAOneScoringFactory taskFactory) {
		super(targetList, library, diaFile, fastaFile, percolatorFiles, taskFactory);
	}
	
	// internal
	public VariantXCorDIAJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, Optional<LibraryInterface> library, File diaFile, StripeFileInterface diaFileReader, File fastaFile, PercolatorExecutionData percolatorFiles,
			XCorDIAOneScoringFactory taskFactory) {
		super(targetList, library, diaFile, diaFileReader, fastaFile, percolatorFiles, taskFactory);
	}

	public VariantXCorDIAJobData updateTaskFactory(XCorDIAOneScoringFactory taskFactory) {
		return new VariantXCorDIAJobData(getTargetList(), getLibrary(), getDiaFile(), getDiaFileReader(), getFastaFile(), getPercolatorFiles(), taskFactory);
	}

	public File getLocalizationFile() {
		String absolutePath = getPrefixFromOutput(getPercolatorFiles().getPeptideOutputFile());
		return new File(absolutePath+".localizations.txt");
	}
	
	@Override
	public void writeToXML(Document doc, Element parentElement) {
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);

		XMLUtils.writeTag(doc, rootElement, "diaFile", getDiaFile().getAbsolutePath());
		XMLUtils.writeTag(doc, rootElement, "fastaFile", getFastaFile().getAbsolutePath());
		if (getLibrary().isPresent()&&getLibrary().get() instanceof LibraryFile) {
			XMLUtils.writeTag(doc, rootElement, "library", ((LibraryFile) getLibrary().get()).getFile().getAbsolutePath());
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

	
	public static VariantXCorDIAJobData readFromXML(Document doc, Element rootElement) {
		if (!rootElement.getTagName().equals(VariantXCorDIAJobData.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+VariantXCorDIAJobData.class.getSimpleName()+"]");
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

		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, readParams);
		BackgroundFrequencyInterface background;
		try {
			background=BackgroundFrequencyCalculator.generateBackground(stripefile);
		} catch (IOException|SQLException|DataFormatException e) {
			throw new EncyclopediaException("Error creating background calculator for VariantXCorDIA", e);
		}
		//BackgroundFrequencyInterface background=new UnitBackgroundFrequencyCalculator(0.01f); //FIXME REMOVE, TESTING ONLY
		VariantXCorDIAOneScoringFactory factory=new VariantXCorDIAOneScoringFactory(readParams, background, new LinkedBlockingQueue<ModificationLocalizationData>());
		return new VariantXCorDIAJobData(Optional.ofNullable(targetList), Optional.ofNullable(libraryObject), diaFile, fastaFile, percolatorData, factory);
	}
}
