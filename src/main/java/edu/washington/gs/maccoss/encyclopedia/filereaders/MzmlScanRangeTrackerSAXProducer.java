package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.xml.parsers.SAXParserFactory;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ScanRangeTracker;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.ProgressInputStream;

public class MzmlScanRangeTrackerSAXProducer extends DefaultHandler implements Runnable {
	private final File mzMLFile;
	private final ScanRangeTracker retentionTimesByStripe=new ScanRangeTracker();
	private final SearchParameters parameters;

	public MzmlScanRangeTrackerSAXProducer(File mzMLFile, SearchParameters parameters) {
		this.mzMLFile=mzMLFile;
		this.parameters=parameters;
	}

	public ScanRangeTracker getRetentionTimesByStripe() {
		return retentionTimesByStripe;
	}

	@Override
	public void run() {
		try {
			final ProgressInputStream stream=new ProgressInputStream(new FileInputStream(mzMLFile));
			final long length=mzMLFile.length();

			stream.addChangeListener(new ChangeListener() {
				int lastUpdate=0;

				@Override
				public void stateChanged(ChangeEvent e) {
					int floor=(int)((stream.getProgress()*100L)/length);
					if (floor>lastUpdate) {
						Logger.logLine("Parsed "+floor+"%");
						lastUpdate=floor;
					}
				}
			});

			SAXParserFactory.newInstance().newSAXParser().parse(stream, this);
		} catch (SAXTerminatorException sax) {
			// ok! this just is an indicator of early termination
			
		} catch (Exception e) {
			Logger.errorLine("Mzml reading failed!");
			Logger.errorException(e);
		}
	}

	private final ArrayList<String> tagList=new ArrayList<String>();
	private Integer msLevel=null;
	private String spectrumRef=null;

	private Float scanStartTime=null;
	private Float isolationWindowTarget=null;
	private Float isolationWindowLowerOffset=null;
	private Float isolationWindowUpperOffset=null;
	private Float scanWindowLowerLimit=null;
	private Float scanWindowUpperLimit=null;
	
	private Float selectedIon=null;

	@Override
	public void startElement(String uri, String localName, String qName, Attributes attributes) throws SAXException {
		if (tagList.size()>0&&"cvParam".equalsIgnoreCase(qName)) {
			if ("spectrum".equalsIgnoreCase(tagList.get(tagList.size()-1))) {
				if ("ms level".equalsIgnoreCase(attributes.getValue("name"))) {
					msLevel=Integer.parseInt(attributes.getValue("value"));
				}
				
			} else if ("scan".equalsIgnoreCase(tagList.get(tagList.size()-1))) {
				if ("scan start time".equalsIgnoreCase(attributes.getValue("name"))) {
					float multiplier;
					String unit=attributes.getValue("unitName");
					if ("second".equalsIgnoreCase(unit)) {
						multiplier=1.0f;
					} else if ("minute".equalsIgnoreCase(unit)) {
						multiplier=60.0f;
					} else if ("hour".equalsIgnoreCase(unit)) {
						multiplier=360.0f;
					} else if ("millisecond".equalsIgnoreCase(unit)) {
						multiplier=0.001f;
					} else {
						throw new EncyclopediaException("Unexpected time unit: "+unit);
					}

					scanStartTime=multiplier*Float.parseFloat(attributes.getValue("value"));
				}

			} else if ("isolationWindow".equalsIgnoreCase(tagList.get(tagList.size()-1))) {
				if ("isolation window target m/z".equalsIgnoreCase(attributes.getValue("name"))) {
					isolationWindowTarget=Float.parseFloat(attributes.getValue("value"));
				} else if ("isolation window lower offset".equalsIgnoreCase(attributes.getValue("name"))) {
					isolationWindowLowerOffset=Float.parseFloat(attributes.getValue("value"));
				} else if ("isolation window upper offset".equalsIgnoreCase(attributes.getValue("name"))) {
					isolationWindowUpperOffset=Float.parseFloat(attributes.getValue("value"));
				}

			} else if ("scanWindow".equalsIgnoreCase(tagList.get(tagList.size()-1))) {
				if ("scan window lower limit".equalsIgnoreCase(attributes.getValue("name"))) {
					scanWindowLowerLimit=Float.parseFloat(attributes.getValue("value"));
				} else if ("scan window upper limit".equalsIgnoreCase(attributes.getValue("name"))) {
					scanWindowUpperLimit=Float.parseFloat(attributes.getValue("value"));
				}

			} else if ("selectedIon".equalsIgnoreCase(tagList.get(tagList.size()-1))) {
				if ("selected ion m/z".equalsIgnoreCase(attributes.getValue("name"))) {
					selectedIon=Float.parseFloat(attributes.getValue("value"));
				}
			}
		} else if ("precursor".equalsIgnoreCase(qName)) {
			spectrumRef=attributes.getValue("spectrumRef");
		}

		tagList.add(qName);
	}

	@Override
	public void endElement(String uri, String localName, String qName) throws SAXException {
		if ("spectrum".equalsIgnoreCase(qName)) {
			if (spectrumRef==null&&msLevel<=1) {
				Range range=new Range(scanWindowLowerLimit, scanWindowUpperLimit);
				retentionTimesByStripe.addPrecursor(range, scanStartTime);
			} else {
				if (spectrumRef==null) spectrumRef="Unknown";

				// Issue 3
				if (isolationWindowTarget==null||isolationWindowLowerOffset==null||isolationWindowUpperOffset==null) {
					if (parameters.getPrecursorWindowSize()>0f&&selectedIon!=null) {
						isolationWindowTarget=selectedIon;
						isolationWindowLowerOffset=parameters.getPrecursorWindowSize()/2.0f;
						isolationWindowUpperOffset=parameters.getPrecursorWindowSize()/2.0f;
					} else {
						Logger.errorLine("Isolation window information missing without precursor window size supplied!");
					}
				}

				Range range=new Range(isolationWindowTarget-isolationWindowLowerOffset, isolationWindowTarget+isolationWindowUpperOffset);
				boolean keepGoing=retentionTimesByStripe.addRange(range, scanStartTime);
				if (!keepGoing) {
					throw new SAXTerminatorException();
				}
			}
			
			spectrumRef=null;

			scanStartTime=null;
			isolationWindowTarget=null;
			isolationWindowLowerOffset=null;
			isolationWindowUpperOffset=null;
			scanWindowLowerLimit=null;
			scanWindowUpperLimit=null;
			
			selectedIon=null;
		}

		tagList.remove(tagList.size()-1);
	}

	@Override
	public void characters(char[] ch, int start, int length) throws SAXException {
	}

	@Override
	public void endDocument() throws SAXException {
	}
	
	public class SAXTerminatorException extends EncyclopediaException {
		public SAXTerminatorException() {
			super("Terminate Early");
		}
	}
}
