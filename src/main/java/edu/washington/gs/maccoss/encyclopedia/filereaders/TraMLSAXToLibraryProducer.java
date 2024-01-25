package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Optional;
import java.util.Map.Entry;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParserFactory;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.io.ProgressInputStream;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class TraMLSAXToLibraryProducer extends DefaultHandler implements Runnable {
	private static final float UNSET_RETENTION_TIME=-Float.MAX_VALUE;
	private static final float UNSET_PEAK_INTENSITY=-Float.MAX_VALUE;
	private final File tramlFile;
	private final SearchParameters parameters;
	private Throwable error;
	ArrayList<LibraryEntry> entries=new ArrayList<>();
	
	public TraMLSAXToLibraryProducer(File tramlFile, SearchParameters parameters) {
		super();
		this.tramlFile=tramlFile;
		this.parameters = parameters;
	}

	@Override
	public void run() {
		try {
			final ProgressInputStream stream = new ProgressInputStream(new FileInputStream(tramlFile));
			final long length = tramlFile.length();
	
			stream.addChangeListener(new ChangeListener() {
				int lastUpdate = 0;
	
				@Override
				public void stateChanged(ChangeEvent e) {
					int floor = (int) ((stream.getProgress() * 100L) / length);
					if (floor > lastUpdate) {
						Logger.logLine("Parsed " + floor + "%");
						lastUpdate = floor;
					}
				}
			});
	
			SAXParserFactory.newInstance().newSAXParser().parse(stream, this);

		} catch (SAXException e) {
			Logger.errorLine("XML Error parsing TraML:");
			Logger.errorException(e);
			this.error = e;

		} catch (ParserConfigurationException e) {
			Logger.errorLine("XML Error parsing TraML:");
			Logger.errorException(e);
			this.error = e;

		} catch (IOException e) {
			Logger.errorLine("IO Error parsing TraML:");
			Logger.errorException(e);
			this.error = e;
		} 
	}

	private final StringBuilder dataSB=new StringBuilder();
	private final ArrayList<String> tagList=new ArrayList<String>();
	
	private String lastPeptideRef=null;
	private String lastPeptideSequence=null;
	private String lastProteinRef=null;
	private byte lastChargeState=0;
	private float lastRetentionTime=UNSET_RETENTION_TIME;
	private ArrayList<ModificationObject> ptms=new ArrayList<>();
	
	private float lastIntensity=UNSET_PEAK_INTENSITY;
	private double lastIsolationWindowTarget=0.0;
	
	HashMap<String, ArrayList<Peak>> transitions=new HashMap<>();
	HashMap<String, PrecursorObject> precursors=new HashMap<>();
	
	@Override
	public void startElement(String uri, String localName, String qName, Attributes attributes) throws SAXException {
		//System.out.println(General.formatCellToWidth("<"+qName+">", tagList.size()*4+qName.length()+3, false));
		dataSB.setLength(0);
		if (tagList.size()>0&&"cvParam".equalsIgnoreCase(qName)) {
			if ("Transition".equalsIgnoreCase(getPreviousElementTag())) {
				if ("product ion intensity".equalsIgnoreCase(attributes.getValue("name"))) {
					lastIntensity=Float.parseFloat(attributes.getValue("value"));
				}
			} else if ("Product".equalsIgnoreCase(getPreviousElementTag())) {
				if ("isolation window target m/z".equalsIgnoreCase(attributes.getValue("name"))) {
					lastIsolationWindowTarget=Double.parseDouble(attributes.getValue("value"));
				}

			} else if ("Peptide".equalsIgnoreCase(getPreviousElementTag())) {
				if ("charge state".equalsIgnoreCase(attributes.getValue("name"))) {
					lastChargeState=Byte.parseByte(attributes.getValue("value"));
				}
				
			} else if ("RetentionTime".equalsIgnoreCase(getPreviousElementTag())) {
				if ("normalized retention time".equalsIgnoreCase(attributes.getValue("name"))) {
					lastRetentionTime=Float.parseFloat(attributes.getValue("value"))*60;
				} else if ("local retention time".equalsIgnoreCase(attributes.getValue("name"))) {
					lastRetentionTime=Float.parseFloat(attributes.getValue("value"))*60;
				} else if ("predicted retention time".equalsIgnoreCase(attributes.getValue("name"))) {
					if (lastRetentionTime==UNSET_RETENTION_TIME) {
						lastRetentionTime=Float.parseFloat(attributes.getValue("value"))*60;
					}
				}
			}
		} else if ("Transition".equalsIgnoreCase(qName)) {
			lastPeptideRef=attributes.getValue("peptideRef");
		} else if ("Peptide".equalsIgnoreCase(qName)) {
			lastPeptideRef=attributes.getValue("id");
			lastPeptideSequence=attributes.getValue("sequence");
		} else if ("ProteinRef".equalsIgnoreCase(qName)) {
			lastProteinRef=attributes.getValue("id");
		} else if ("Modification".equalsIgnoreCase(qName)) {
			int location=Integer.parseInt(attributes.getValue("location"));
			double mass=Double.parseDouble(attributes.getValue("monoisotopicMassDelta"));
			ptms.add(new ModificationObject(location, mass));
		}
		tagList.add(qName);
	}

	@Override
	public void endElement(String uri, String localName, String qName) throws SAXException {
		//System.out.println(General.formatCellToWidth("</"+qName+">", (tagList.size()-1)*4+qName.length()+4, false));
		if ("Transition".equalsIgnoreCase(qName)) {
			if (lastPeptideRef!=null) {
				ArrayList<Peak> peaks=transitions.get(lastPeptideRef);
				if (peaks==null) {
					peaks=new ArrayList<>();
					transitions.put(lastPeptideRef, peaks);
				}
				if (lastIntensity==UNSET_PEAK_INTENSITY) {
					// transition lists don't necessarily have intensities, so this forces there to be something tiny
					lastIntensity=Float.MIN_VALUE;
				}
				peaks.add(new Peak(lastIsolationWindowTarget, lastIntensity));
			}
			lastPeptideRef=null;
			lastIsolationWindowTarget=0.0;
			lastIntensity=UNSET_PEAK_INTENSITY;
		} else if ("Peptide".equalsIgnoreCase(qName)) {
			if (lastPeptideRef!=null) {
				String peptideModSeq=getPeptideModSeq(lastPeptideSequence, ptms);
				PrecursorObject precursor=new PrecursorObject(peptideModSeq, lastProteinRef, lastChargeState, lastRetentionTime);
				precursors.put(lastPeptideRef, precursor);				
			}
			lastPeptideRef=null;
			lastPeptideSequence=null;
			lastProteinRef=null;
			lastChargeState=0;
			lastRetentionTime=UNSET_RETENTION_TIME;
			ptms.clear();
		}
		tagList.remove(tagList.size()-1);
	}

	@Override
	public void characters(char[] ch, int start, int length) throws SAXException {
		// THIS PARSER DOESN'T USE DATA
		//dataSB.append(ch, start, length);
	}

	@Override
	public void endDocument() throws SAXException {
		String sourceFile = tramlFile.getName();
		AminoAcidConstants aaConstants=parameters.getAAConstants();
		
		for (Entry<String, PrecursorObject> entry : precursors.entrySet()) {
			String peptideRef=entry.getKey();
			PrecursorObject precursor=entry.getValue();
			ArrayList<Peak> transitionList=transitions.get(peptideRef);
			if (transitionList==null) {
				Logger.errorLine("Found precusor ("+precursor.peptideModSeq+") but no peaks! Skipping entry.");
				continue;
			}

			Triplet<double[], float[], Optional<float[]>> peakArrays=Peak.toArrays(transitionList);
			
			double precursorMZ=aaConstants.getChargedMass(precursor.peptideModSeq, precursor.chargeState);

			HashSet<String> accessions=new HashSet<>();
			accessions.add(precursor.proteinReference);
			LibraryEntry libEntry=new LibraryEntry(sourceFile, accessions, precursorMZ, precursor.chargeState, precursor.peptideModSeq, 1, precursor.retentionTime, 0.0f, peakArrays.x, peakArrays.y, Optional.empty(), // FIXME add ion mobility to parser
					aaConstants);
			entries.add(libEntry);
		}
	}
	
	private static class PrecursorObject {
		private final String peptideModSeq;
		private final String proteinReference;
		private final byte chargeState;
		private final float retentionTime;
		public PrecursorObject(String peptideModSeq, String proteinReference, byte chargeState, float retentionTime) {
			this.peptideModSeq = peptideModSeq;
			this.proteinReference = proteinReference;
			this.chargeState = chargeState;
			this.retentionTime = retentionTime;
		}
	}
	
	private static class ModificationObject {
		private final int location;
		private final double mass;
		public ModificationObject(int location, double mass) {
			this.location = location;
			this.mass = mass;
		}
	}

	private String getPreviousElementTag() {
		return tagList.get(tagList.size() - 1);
	}
	
	private static String getPeptideModSeq(String sequence, ArrayList<ModificationObject> ptms) {
		String[] aas=new String[sequence.length()];
		double[] ptmMasses=new double[aas.length];
		for (int i=0; i<aas.length; i++) {
			aas[i]=Character.toString(sequence.charAt(i));
		}
		for (ModificationObject ptm : ptms) {
			Double modMass=ptm.mass;
			int index=ptm.location-1;
			if (index==aas.length) index=aas.length-1;
			if (index<0) index=0;
			ptmMasses[index]+=modMass;
		}
		for (int i = 0; i < aas.length; i++) {
			if (ptmMasses[i]!=0) {
				aas[i]=aas[i]+"["+(ptmMasses[i]>=0?"+":"")+ptmMasses[i]+"]";
			}
		}
		StringBuilder sb=new StringBuilder();
		for (int i=0; i<aas.length; i++) {
			sb.append(aas[i]);
		}
		String peptideModSeq=sb.toString();
		return peptideModSeq;
	}


	public boolean hadError() {
		return null != error;
	}

	public Throwable getError() {
		return error;
	}
	
	public ArrayList<LibraryEntry> getEntries() {
		return entries;
	}
}
