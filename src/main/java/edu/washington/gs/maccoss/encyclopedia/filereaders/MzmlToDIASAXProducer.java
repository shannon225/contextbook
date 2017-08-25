package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Base64;
import java.util.HashMap;
import java.util.concurrent.BlockingQueue;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.xml.parsers.SAXParserFactory;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.ProgressInputStream;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;

public class MzmlToDIASAXProducer extends DefaultHandler implements Runnable {
	private final File mzMLFile;
	private final BlockingQueue<MzmlBlock> mzmlBlockQueue;
	private final SearchParameters parameters;
	private final HashMap<Range, TFloatArrayList> retentionTimesByStripe=new HashMap<Range, TFloatArrayList>();

	private Throwable error;

	public MzmlToDIASAXProducer(File mzMLFile, BlockingQueue<MzmlBlock> mzmlBlockQueue, SearchParameters parameters) {
		this.mzMLFile=mzMLFile;
		this.mzmlBlockQueue=mzmlBlockQueue;
		this.parameters=parameters;
	}

	public HashMap<Range, TFloatArrayList> getRetentionTimesByStripe() {
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

			// Just for safety, ensure that if we get this far the consumer(s) of the queue will finish
			// If parse() already put a block, this will never be consumed, but it can't hurt
			mzmlBlockQueue.put(MzmlBlock.POISON_BLOCK);
		} catch (Throwable t) {
			Logger.errorLine("Mzml reading failed!");
			Logger.errorException(t);

			this.error = t;
		}
	}

	public boolean hadError() {
		return null != error;
	}

	public Throwable getError() {
		return error;
	}

	public String getMzMLID() {
		return mzML_ID;
	}

	private final ArrayList<PrecursorScan> precursors=new ArrayList<PrecursorScan>();
	private final ArrayList<Stripe> stripes=new ArrayList<Stripe>();

	private final ArrayList<String> tagList=new ArrayList<String>();
	private String mzML_ID=null;
	private String spectrumName=null;
	private Integer spectrumIndex=null;
	private Integer msLevel=null;
	private String spectrumRef=null;

	private Float scanStartTime=null;
	private Float isolationWindowTarget=null;
	private Float isolationWindowLowerOffset=null;
	private Float isolationWindowUpperOffset=null;

	private boolean compress=false;
	private Precision encoding=null;
	private Boolean isMZ=null;
	private double[] massArray=null;
	private float[] intensityArray=null;
	private final StringBuilder dataSB=new StringBuilder();
	
	private Float selectedIon=null;
	private Byte selectedCharge=null;

	private boolean isSkipSpectrumWithBadEncoding = false;
	private int numSkippedWithBadEncoding = 0;
	private static final int MAX_BAD_ENCODING_SKIPPED = 1000;
	
	@Override
	public void startElement(String uri, String localName, String qName, Attributes attributes) throws SAXException {
		dataSB.setLength(0);
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

			} else if (tagList.size()>2&&"binaryDataArray".equals(tagList.get(tagList.size()-1))) {
				if ("binaryDataArrayList".equals(tagList.get(tagList.size()-2))) {
					if ("spectrum".equals(tagList.get(tagList.size()-3))) {
						if ("64-bit float".equalsIgnoreCase(attributes.getValue("name"))) {
							encoding=Precision.FLOAT64BIT;
						} else if ("64-bit integer".equalsIgnoreCase(attributes.getValue("name"))) {
							encoding=Precision.INT64BIT;
						} else if ("32-bit float".equalsIgnoreCase(attributes.getValue("name"))) {
							encoding=Precision.FLOAT32BIT;
						} else if ("32-bit integer".equalsIgnoreCase(attributes.getValue("name"))) {
							encoding=Precision.INT32BIT;
						} else if ("zlib compression".equalsIgnoreCase(attributes.getValue("name"))) {
							compress=true;
						} else if ("m/z array".equalsIgnoreCase(attributes.getValue("name"))) {
							isMZ=true;
						} else if ("intensity array".equalsIgnoreCase(attributes.getValue("name"))) {
							isMZ=false;
						}
					}
				}
				
			} else if ("selectedIon".equalsIgnoreCase(tagList.get(tagList.size()-1))) {
				if ("selected ion m/z".equalsIgnoreCase(attributes.getValue("name"))) {
					selectedIon=Float.parseFloat(attributes.getValue("value"));
				} else if ("charge state".equalsIgnoreCase(attributes.getValue("name"))) {
					selectedCharge=Byte.parseByte(attributes.getValue("value"));
				}
			}
		} else if ("precursor".equalsIgnoreCase(qName)) {
			spectrumRef=attributes.getValue("spectrumRef");
		} else if ("mzML".equalsIgnoreCase(qName)) {
			mzML_ID=attributes.getValue("id");
		} else if ("spectrum".equalsIgnoreCase(qName)) {
			spectrumName=attributes.getValue("id");
			spectrumIndex=Integer.parseInt(attributes.getValue("index"));
		}

		tagList.add(qName);
	}

	@Override
	public void endElement(String uri, String localName, String qName) throws SAXException {
		if ("spectrum".equalsIgnoreCase(qName)) {
			
			if (isSkipSpectrumWithBadEncoding){
				
				Logger.errorLine("Skipping spectrum #" + spectrumIndex + ", '" + spectrumName + "': bad binary data encoding.");
				
				//reset variable fields
				spectrumName=null;
				spectrumIndex=null;
				spectrumRef=null;

				scanStartTime=null;
				isolationWindowTarget=null;
				isolationWindowLowerOffset=null;
				isolationWindowUpperOffset=null;

				compress=false;
				encoding=null;
				isMZ=null;
				dataSB.setLength(0);
				massArray=null;
				intensityArray=null;
				
				selectedIon=null;
				
				//reset skipping flag and increment counter
				isSkipSpectrumWithBadEncoding = false;
				numSkippedWithBadEncoding++;
				
				return;
			}
			
			if (spectrumRef==null&&msLevel<=1) {
				if (parameters.getPrecursorOffsetPPM()!=0.0) {
					double[] deltaArray=General.multiply(massArray, parameters.getPrecursorOffsetPPM()/1000000.0);
					massArray=General.subtract(massArray, deltaArray);
				}
				precursors.add(new PrecursorScan(spectrumName, spectrumIndex, scanStartTime, massArray, intensityArray));

			} else {
				if (spectrumRef==null) spectrumRef="Unknown";
				
				if (parameters.getFragmentOffsetPPM()!=0.0) {
					double[] deltaArray=General.multiply(massArray, parameters.getFragmentOffsetPPM()/1000000.0);
					massArray=General.subtract(massArray, deltaArray);
				}
				
				//Issue 3
				if (isolationWindowTarget==null||isolationWindowLowerOffset==null||isolationWindowUpperOffset==null) {
					if (parameters.getPrecursorWindowSize()>0f&&selectedIon!=null) {
						isolationWindowTarget=selectedIon;
						isolationWindowLowerOffset=parameters.getPrecursorWindowSize()/2.0f;
						isolationWindowUpperOffset=parameters.getPrecursorWindowSize()/2.0f;
					} else {
						Logger.errorLine("Isolation window information missing without precursor window size supplied!");
					}
				}
				
				byte charge=0;
				if (selectedCharge!=null) {
					charge=selectedCharge;
				}
				
				Stripe stripe=new Stripe(spectrumName, spectrumRef, spectrumIndex, scanStartTime, isolationWindowTarget-isolationWindowLowerOffset+(float)parameters.getPrecursorIsolationMargin(), isolationWindowTarget+isolationWindowUpperOffset-(float)parameters.getPrecursorIsolationMargin(),
						massArray, intensityArray, charge);
				stripes.add(stripe);
				Range range=stripe.getRange();
				TFloatArrayList stripeRTs=retentionTimesByStripe.get(range);
				if (stripeRTs==null) {
					stripeRTs=new TFloatArrayList();
					retentionTimesByStripe.put(range, stripeRTs);
				}
				stripeRTs.add(scanStartTime);
				
			}

			try {
				if (precursors.size()>100||stripes.size()>1000) {
					mzmlBlockQueue.put(new MzmlBlock(precursors, stripes));
					precursors.clear();
					stripes.clear();
				}
			} catch (InterruptedException ie) {
				Logger.errorLine("Mzml reading interrupted!");
				Logger.errorException(ie);
			}

			spectrumName=null;
			spectrumIndex=null;
			spectrumRef=null;

			scanStartTime=null;
			isolationWindowTarget=null;
			isolationWindowLowerOffset=null;
			isolationWindowUpperOffset=null;

			compress=false;
			encoding=null;
			isMZ=null;
			dataSB.setLength(0);
			massArray=null;
			intensityArray=null;
			
			selectedIon=null;
		}

		if (tagList.size()>3&&"binary".equals(tagList.get(tagList.size()-1))) {
			if ("binaryDataArray".equals(tagList.get(tagList.size()-2))) {
				if ("binaryDataArrayList".equals(tagList.get(tagList.size()-3))) {
					if ("spectrum".equals(tagList.get(tagList.size()-4))) {
						try {
							// isMZ required to be set every time!
							// data can be null (0 length array)
							if (isMZ!=null) {
								if (encoding==null) encoding=Precision.FLOAT64BIT;
								Number[] values=getBinaryDataAsNumberArray(dataSB.toString(), compress, encoding);
								if (isMZ) {
									massArray=ByteConverter.toDoubleArray(values);
								} else {
									intensityArray=ByteConverter.toFloatArray(values);
								}

								compress=false;
								encoding=null;
								dataSB.setLength(0);
								isMZ=null;
							}
						} catch (IOException e) {
							throw new EncyclopediaException("Error parsing binary data from "+General.toString(tagList), e);
						}
					}
				}
			}
		}

		tagList.remove(tagList.size()-1);
	}

	@Override
	public void characters(char[] ch, int start, int length) throws SAXException {
		if ("binary".equals(tagList.get(tagList.size()-1))) {
			dataSB.append(ch, start, length);
		}
	}

	@Override
	public void endDocument() throws SAXException {
		try {
			mzmlBlockQueue.put(new MzmlBlock(precursors, stripes));
			mzmlBlockQueue.put(MzmlBlock.POISON_BLOCK);
			if (numSkippedWithBadEncoding > MAX_BAD_ENCODING_SKIPPED) {
				throw new EncyclopediaException("Skipped too many spectra as a result of bad binary data encoding! File is invalid.");
			} else if (numSkippedWithBadEncoding > 0){
				Logger.errorLine("Skipped " + numSkippedWithBadEncoding + " spectra because of bad binary data encoding.");
			}
		} catch (InterruptedException ie) {
			Logger.errorLine("Mzml reading interrupted!");
			Logger.errorException(ie);
		}
	}

	public Number[] getBinaryDataAsNumberArray(String binaryString, boolean decompress, Precision precision) throws IOException {
		if (binaryString==null||binaryString.length()==0) return new Number[0];
		
		byte[] binArray = null;
		try {
			binArray=Base64.getDecoder().decode(binaryString);
		} catch (IllegalArgumentException e) {
			isSkipSpectrumWithBadEncoding = true;
			//System.out.println(binaryString);
			//throw new IllegalArgumentException(e);
			return new Number[0];	//dummy output
		}
		
		ByteBuffer bbuf=ByteBuffer.allocate(binArray.length);
		bbuf.put(binArray);
		binArray=bbuf.order(ByteOrder.LITTLE_ENDIAN).array();

		byte[] data;
		if (decompress) {
			try {
				data=CompressionUtils.decompress(binArray);
			} catch (IllegalStateException e) {
				//throw new EncyclopediaException("Error parsing binary data from "+binaryString, e);
				isSkipSpectrumWithBadEncoding = true;
				return new Number[0];	//dummy output
			}
		} else {
			data=binArray;
		}

		return convertData(data, precision);
	}

	private Number[] convertData(byte[] data, Precision prec) {
		int step;
		switch (prec) {
			case FLOAT64BIT:
			case INT64BIT:
				step=8;
				break;
			case FLOAT32BIT:
			case INT32BIT:
				step=4;
				break;
			default:
				step=-1;
		}
		Number[] resultArray=new Number[data.length/step];
		ByteBuffer bb=ByteBuffer.wrap(data);
		bb.order(ByteOrder.LITTLE_ENDIAN); // the order in mzML is LITTLE_ENDIAN
		for (int indexOut=0; indexOut<data.length; indexOut+=step) {
			Number num;
			switch (prec) {
				case FLOAT64BIT:
					num=bb.getDouble(indexOut);
					break;
				case INT64BIT:
					num=bb.getLong(indexOut);
					break;
				case FLOAT32BIT:
					num=bb.getFloat(indexOut);
					break;
				case INT32BIT:
					num=bb.getInt(indexOut);
					break;
				default:
					num=null;
			}
			resultArray[indexOut/step]=num;
		}
		return resultArray;
	}

	public enum Precision {
		FLOAT32BIT, FLOAT64BIT, INT32BIT, INT64BIT
	}
}
