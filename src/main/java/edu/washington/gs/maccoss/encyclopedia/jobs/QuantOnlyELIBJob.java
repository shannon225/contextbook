package edu.washington.gs.maccoss.encyclopedia.jobs;

import edu.washington.gs.maccoss.encyclopedia.SearchToBLIB;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLObject;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import java.io.File;
import java.util.ArrayList;

/**
 * {@code QuantOnlyJob} runs the quant-only workflow on a subset of the
 * jobs with results recorded in an alignment-only library ("ALIB").
 * When creating a new job, it can either use an existing alignment-only
 * file, or if not specified, the last {@link AlignmentOnlyLibraryJob}
 * from the provided processor.
 *
 * This job will always process all the {@code SearchJob} objects in the
 * provided processor's queue *at the time it runs*; this may severely
 * limit the utility of this class.
 */
public class QuantOnlyELIBJob implements WorkerJob, XMLObject {
	protected final File destFile, alignmentFile;

	protected final JobProcessor processor;

	public QuantOnlyELIBJob(File destFile, File alignmentFile, JobProcessor processor) {
		this.destFile = destFile;
		this.processor=processor;

		if (null != alignmentFile) {
			this.alignmentFile = alignmentFile;
		} else {
			final AlignmentOnlyLibraryJob alignmentJob = lastAlignmentJob(processor);
			if (null == alignmentJob) {
				throw new IllegalArgumentException("No alignment-only job was found in the queue! Did you mean to specify an existing alignment-only file?");
			}
			this.alignmentFile = alignmentJob.destFile;
		}
	}

	private AlignmentOnlyLibraryJob lastAlignmentJob(JobProcessor processor) {
		final ArrayList<WorkerJob> queue = processor.getQueue();

		for (int i = queue.size(); --i >= 0; ) {
			final WorkerJob job = queue.get(i);
			if (job instanceof AlignmentOnlyLibraryJob) {
				return (AlignmentOnlyLibraryJob) job;
			}
		}

		return null;
	}

	@Override
	public String getJobTitle() {
		return "Quantify peptides using " + alignmentFile.getAbsolutePath() + " (writes to " + destFile.getAbsolutePath() + ")";
	}

	@Override
	public void runJob(ProgressIndicator progress) throws Exception {
		ArrayList<SearchJobData> jobData=new ArrayList<SearchJobData>();
		for (WorkerJob job : processor.getQueue()) {
			if (job instanceof SearchJob) {
				jobData.add(((SearchJob)job).getSearchData());
			}
		}

		Logger.logLine("Found "+jobData.size()+" jobs in the queue to combine...");

		if (!jobData.isEmpty()) {
			final SearchJobData representative = jobData.iterator().next();
			Logger.logLine("Extracting representative search parameters from job " + representative.getOriginalDiaFileName());

			SearchToBLIB.convertElibQuantOnly(
					progress,
					jobData,
					destFile,
					alignmentFile,
					representative.getParameters()
			);
		}
	}

	@Override
	public void writeToXML(Document doc, Element parentElement) {
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);

		XMLUtils.writeTag(doc, rootElement, "elibFile", destFile.getAbsolutePath());
		XMLUtils.writeTag(doc, rootElement, "alignmentFrom", alignmentFile.getAbsolutePath());
	}

	public static QuantOnlyELIBJob readFromXML(Document doc, Element rootElement, JobProcessor processor) {
		if (!rootElement.getTagName().equals(QuantOnlyELIBJob.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+QuantOnlyELIBJob.class.getSimpleName()+"]");
		}
		File elibFile=null;
		File alignmentFile=null;

		NodeList nodes=rootElement.getChildNodes();

		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
			if (node.getNodeType() == Node.ELEMENT_NODE) {
				Element element = (Element) node;
				if ("elibFile".equals(element.getTagName())) {
					elibFile=new File(element.getTextContent());
				} else if ("alignmentFrom".equals(element.getTagName())) {
					alignmentFile = new File(element.getTextContent());
				}
			}
		}

		if (elibFile==null) throw new EncyclopediaException("Found null elibFile in "+rootElement.getTagName());

		// It's OK if the alignment file is null here; the constructor will check for a previous alignment-only
		// job in the processor's queue.
		return new QuantOnlyELIBJob(elibFile, alignmentFile, processor);
	}
}
