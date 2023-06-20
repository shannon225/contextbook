package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashSet;
import java.util.Iterator;
import java.util.stream.Collectors;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import com.google.common.collect.ImmutableSortedSet;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLObject;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import gnu.trove.map.hash.TIntIntHashMap;

public class FastaPeptideEntry implements Comparable<FastaPeptideEntry>, XMLObject {
	private final String filename;
	private final HashSet<String> accessions=new HashSet<String>(); // not immutable!
	private final String sequence;
	private final String sequenceWithModsStripped;

	public FastaPeptideEntry(String filename, HashSet<String> accessions, String sequence) {
		this.filename=filename;
		this.accessions.addAll(accessions);
		this.sequence=sequence;
		this.sequenceWithModsStripped=stripMods(sequence);
	}

	public FastaPeptideEntry(String filename, String accession, String sequence) {
		this.filename=filename;
		accessions.add(accession);
		this.sequence=sequence;
		this.sequenceWithModsStripped=stripMods(sequence);
	}

	public FastaPeptideEntry(String sequence) {
		this("Unknown File", "Unknown Annotation", sequence);
	}
	
	public FastaEntryInterface toFastaEntryInterface() {
		return new FastaEntry(filename, PSMData.accessionsToString(accessions), sequence);
	}

	@Override
	public void writeToXML(Document doc, Element parentElement) {
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);

		XMLUtils.writeTag(doc, rootElement, "filename", filename);
		XMLUtils.writeTag(doc, rootElement, "sequence", sequence);
		XMLUtils.writeTag(doc, rootElement, "accessions", PSMData.accessionsToString(accessions));
	}
	
	public static FastaPeptideEntry readFromXML(Document doc, Element rootElement) {
		if (!rootElement.getTagName().equals(FastaPeptideEntry.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+FastaPeptideEntry.class.getSimpleName()+"]");
		}
		String filename=null;
		String sequence=null;
		String accessions=null;
		
		NodeList nodes=rootElement.getChildNodes();
		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
            if (node.getNodeType() == Node.ELEMENT_NODE) {
                Element element = (Element) node;
                if ("filename".equals(element.getTagName())) {
                	filename=element.getTextContent();
                } else if ("sequence".equals(element.getTagName())) {
                	sequence=element.getTextContent();
                } else if ("accessions".equals(element.getTagName())) {
                	accessions=element.getTextContent();
                }
            }
		}
		
		if (filename==null) throw new EncyclopediaException("Found null filename in "+rootElement.getTagName());
		if (sequence==null) throw new EncyclopediaException("Found null sequence in "+rootElement.getTagName());
		if (accessions==null) throw new EncyclopediaException("Found null accessions in "+rootElement.getTagName());
		
		return new FastaPeptideEntry(filename, PSMData.stringToAccessions(accessions), sequence);
	}
	
	private String stripMods(String sequence) {
		StringBuilder sb=new StringBuilder();
		
		boolean inmod=false;
		for (int i = 0; i < sequence.length(); i++) {
			if (sequence.charAt(i)=='[') {
				inmod=true;
			} else if (sequence.charAt(i)==']') {
				inmod=false;
			} else if (!inmod) {
				sb.append(sequence.charAt(i));
			}
		}
		
		return sb.toString();
	}
	
	public void addAccessions(HashSet<String> accessions) {
		this.accessions.addAll(accessions);
	}

	@Override
	public int compareTo(FastaPeptideEntry o) {
		if (o==null) return 1;

		int c;

		final Iterator<String>
				iter = ImmutableSortedSet.copyOf(getAccessions()).iterator(),
				oIter = ImmutableSortedSet.copyOf(o.getAccessions()).iterator();
		while (iter.hasNext() && oIter.hasNext()) {
			String acc = iter.next(), oAcc = oIter.next();

			c = acc.compareTo(oAcc);
			if (c!=0) return c;
		}

		c = getAccessions().size() - o.getAccessions().size();
		if (c!=0) return c;

		c=filename.compareTo(o.getFilename());
		if (c!=0) return c;
		c=sequence.compareTo(o.getSequence());
		
		return c;
		
	}

	public HashSet<String> getAccessions() {
		return accessions;
	}


	/**
	 * Returns this entry's set of accessions with the given flag
	 * prepended to each. Useful for generating decoys.
	 */
	public HashSet<String> getFlaggedAccessions(String prefix) {
		return getAccessions().stream()
				.map(s -> prefix + s)
				.collect(Collectors.toCollection(HashSet::new));
	}

	public String getFilename() {
		return filename;
	}

	public String getSequence() {
		return sequence;
	}
	
	public String getSequenceWithModsStripped() {
		return sequenceWithModsStripped;
	}

	public void addStatistics(TIntIntHashMap map) {
		FastaEntry.getStatistics(sequence, map);
	}

}
