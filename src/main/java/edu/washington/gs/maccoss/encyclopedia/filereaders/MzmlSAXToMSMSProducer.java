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

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMultimap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.mzml.InstrumentComponent;
import edu.washington.gs.maccoss.encyclopedia.filereaders.mzml.InstrumentId;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.io.ProgressInputStream;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;

public class MzmlSAXToMSMSProducer extends DefaultHandler implements MSMSProducer {
	private final File mzMLFile;
	private final BlockingQueue<MSMSBlock> mzmlBlockQueue;
	private final SearchParameters parameters;
	private final HashMap<Range, TFloatArrayList> retentionTimesByStripe=new HashMap<Range, TFloatArrayList>();
	private final HashMap<Range, TFloatArrayList> ionInjectionTimesByStripe=new HashMap<Range, TFloatArrayList>();
	private final ImmutableMultimap.Builder<String, String> softwareAccessionIdToVersionBuilder = ImmutableMultimap.builder();

	private final ImmutableMultimap.Builder<InstrumentId, InstrumentComponent> instrumentIdToInstrumentComponentBuilder = ImmutableMultimap.builder();
	private ImmutableList.Builder<InstrumentComponent> instrumentComponentsBuilder = ImmutableList.builder();

	private Throwable error;
	private int fraction;

	public MzmlSAXToMSMSProducer(File mzMLFile, int fraction, BlockingQueue<MSMSBlock> mzmlBlockQueue, SearchParameters parameters) {
		this.mzMLFile=mzMLFile;
		this.mzmlBlockQueue=mzmlBlockQueue;
		this.parameters=parameters;
		this.fraction=fraction;
	}

	public HashMap<Range, TFloatArrayList> getRetentionTimesByStripe() {
		return retentionTimesByStripe;
	}
	public HashMap<Range, TFloatArrayList> getIonInjectionTimesByStripe() {
		return ionInjectionTimesByStripe;
	}

	@Override
	public void run() {
		try {
			final ProgressInputStream stream = new ProgressInputStream(new FileInputStream(mzMLFile));
			final long length = mzMLFile.length();

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

			// Just for safety, ensure that if we get this far the consumer(s) of the queue will finish
			// If parse() already put a block, this will never be consumed, but it can't hurt
			putBlock(MSMSBlock.POISON_BLOCK);
		} catch (Throwable t) {
			Logger.errorLine("mzML reading failed!");
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
	private final ArrayList<FragmentScan> stripes=new ArrayList<FragmentScan>();

	private final ArrayList<String> tagList=new ArrayList<String>();
	private String mzML_ID=null;
	private String spectrumName=null;
	private Integer spectrumIndex=null;
	private Integer msLevel=null;
	private String spectrumRef=null;
	private String id=null;
	private HashMap<String, HashMap<String, String>> referenceableParamGroups=new HashMap<>();

	private Float scanStartTime=null;
	private Float ionInjectTime=null;
	private Double isolationWindowTarget=null;
	private Double isolationWindowLowerOffset=null;
	private Double isolationWindowUpperOffset=null;

	private Double scanWindowLowerLimit=null;
	private Double scanWindowUpperLimit=null;

	private boolean compress=false;
	private Precision encoding=null;
	private Boolean isMZ=null;
	private double[] massArray=null;
	private float[] intensityArray=null;
	private Float tic=null;
	private final StringBuilder dataSB=new StringBuilder();
	
	private Double selectedIon=null;
	private Byte selectedCharge=null;

	private boolean isSkipSpectrumWithBadEncoding = false;
	private int numSkippedWithBadEncoding = 0;
	private static final int MAX_BAD_ENCODING_SKIPPED = 1000;

	/**
	 * Field is only used to temporarily store software version when reading.
	 * Use {@link #getSoftwareAccessionIdToVersion()}
 	 */
	private String softwareVersion;

	private InstrumentId.Builder currentInstrumentConfigurationIdBuilder = InstrumentId.builder();
	private InstrumentComponent.Builder currentInstrumentComponentBuilder = InstrumentComponent.builder();

	@Override
	public void startElement(String uri, String localName, String qName, Attributes attributes) throws SAXException {
		dataSB.setLength(0);
		if (tagList.size()>0&&"referenceableParamGroupRef".equalsIgnoreCase(qName)) {
			if ("spectrum".equalsIgnoreCase(tagList.get(tagList.size()-1))) {
				HashMap<String,String> mapping=referenceableParamGroups.get(attributes.getValue("ref"));
				processSpectrumParams(mapping);
			}
			
		} else if (tagList.size()>0&&"cvParam".equalsIgnoreCase(qName)) {
			if ("spectrum".equalsIgnoreCase(tagList.get(tagList.size()-1))) {
				processSpectrumParams(attributes);

			} else if ("referenceableParamGroup".equalsIgnoreCase(tagList.get(tagList.size()-1))) {
				HashMap<String,String> map=referenceableParamGroups.get(id);
				map.put(attributes.getValue("name"), attributes.getValue("value"));
				
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
					} else if ("microsecond".equalsIgnoreCase(unit)) {
						multiplier=0.000001f;
					} else {
						throw new EncyclopediaException("Unexpected time unit: "+unit);
					}

					scanStartTime=multiplier*Float.parseFloat(attributes.getValue("value"));
				} else if ("ion injection time".equalsIgnoreCase(attributes.getValue("name"))) {
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
					} else if ("microsecond".equalsIgnoreCase(unit)) {
						multiplier=0.000001f;
					} else {
						throw new EncyclopediaException("Unexpected time unit: "+unit);
					}

					ionInjectTime=multiplier*Float.parseFloat(attributes.getValue("value")); 
				}
				
			} else if ("isolationWindow".equalsIgnoreCase(tagList.get(tagList.size()-1))) {
				if ("isolation window target m/z".equalsIgnoreCase(attributes.getValue("name"))) {
					isolationWindowTarget=Double.parseDouble(attributes.getValue("value"));
				} else if ("isolation window lower offset".equalsIgnoreCase(attributes.getValue("name"))) {
					isolationWindowLowerOffset=Double.parseDouble(attributes.getValue("value"));
				} else if ("isolation window upper offset".equalsIgnoreCase(attributes.getValue("name"))) {
					isolationWindowUpperOffset=Double.parseDouble(attributes.getValue("value"));
				}
				
			} else if ("scanWindow".equalsIgnoreCase(tagList.get(tagList.size()-1))) {
				if ("scan window lower limit".equalsIgnoreCase(attributes.getValue("name"))) {
					scanWindowLowerLimit=Double.parseDouble(attributes.getValue("value"));
				} else if ("scan window upper limit".equalsIgnoreCase(attributes.getValue("name"))) {
					scanWindowUpperLimit=Double.parseDouble(attributes.getValue("value"));
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
					selectedIon=Double.parseDouble(attributes.getValue("value"));
				} else if ("charge state".equalsIgnoreCase(attributes.getValue("name"))) {
					selectedCharge=Byte.parseByte(attributes.getValue("value"));
				}
			} else if ("software".equalsIgnoreCase(getPreviousElementTag())) {
				softwareAccessionIdToVersionBuilder.put(attributes.getValue("accession"), softwareVersion);
			} else if ("instrumentConfiguration".equalsIgnoreCase(getPreviousElementTag())) {
				currentInstrumentConfigurationIdBuilder
						.setAccession(attributes.getValue("accession"))
						.setName(attributes.getValue("name"));
			} else if (InstrumentComponent.Type.getTypeByName(getPreviousElementTag()).isPresent()) {
				currentInstrumentComponentBuilder
						.setType(InstrumentComponent.Type.getTypeByName(getPreviousElementTag()).get())
						.setCvRef(attributes.getValue("cvRef"))
						.setAccessionId(attributes.getValue("accession"))
						.setName(attributes.getValue("name"));
			}
		} else if ("referenceableParamGroup".equalsIgnoreCase(qName)) {
			id=attributes.getValue("id");
			referenceableParamGroups.put(id, new HashMap<>());
		} else if ("precursor".equalsIgnoreCase(qName)) {
			spectrumRef=attributes.getValue("spectrumRef");
		} else if ("mzML".equalsIgnoreCase(qName)) {
			mzML_ID=attributes.getValue("id");
		} else if ("spectrum".equalsIgnoreCase(qName)) {
			spectrumName=attributes.getValue("id");
			spectrumIndex=Integer.parseInt(attributes.getValue("index"));
		} else if ("software".equalsIgnoreCase(qName)) {
			softwareVersion = attributes.getValue("version");
		} else if (InstrumentComponent.Type.getTypeByName(qName).isPresent()) {
			currentInstrumentComponentBuilder.setOrder(Integer.parseInt(attributes.getValue("order")));
		} else if ("instrumentConfiguration".equalsIgnoreCase(qName)) {
			currentInstrumentConfigurationIdBuilder.setInstrumentConfigurationId(attributes.getValue("id"));
		}

		tagList.add(qName);
	}

	private void processSpectrumParams(Attributes attributes) {
		if ("ms level".equalsIgnoreCase(attributes.getValue("name"))) {
			msLevel=Integer.parseInt(attributes.getValue("value"));
		} else if ("total ion current".equalsIgnoreCase(attributes.getValue("name"))) {
			tic=Float.parseFloat(attributes.getValue("value"));
		}
	}

	private void processSpectrumParams(HashMap<String,String> attributes) {
		if (attributes.containsKey("ms level")) msLevel=Integer.parseInt(attributes.get("ms level"));
		if (attributes.containsKey("total ion current")) tic=Float.parseFloat(attributes.get("total ion current"));
	}

	private String getPreviousElementTag() {
		return tagList.get(tagList.size() - 1);
	}

	/**
	 * @return A mapping of PSI-MS controlled software accession id to version(s).
	 */
	public ImmutableMultimap<String, String> getSoftwareAccessionIdToVersion() {
		return softwareAccessionIdToVersionBuilder.build();
	}

	public ImmutableMultimap<InstrumentId, InstrumentComponent> getInstrumentConfigurations() {
		return instrumentIdToInstrumentComponentBuilder.build();
	}

	@Override
	public void endElement(String uri, String localName, String qName) throws SAXException {
		if ("spectrum".equalsIgnoreCase(qName)) {
			
			if (isSkipSpectrumWithBadEncoding) {
				Logger.errorLine("Skipping spectrum #" + spectrumIndex + ", '" + spectrumName + "': bad binary data encoding.");
				
				//reset variable fields
				spectrumName=null;
				spectrumIndex=null;
				spectrumRef=null;

				scanStartTime=null;
				ionInjectTime=null;
				isolationWindowTarget=null;
				isolationWindowLowerOffset=null;
				isolationWindowUpperOffset=null;
				scanWindowLowerLimit=null;
				scanWindowUpperLimit=null; 

				compress=false;
				encoding=null;
				isMZ=null;
				dataSB.setLength(0);
				massArray=null;
				intensityArray=null;
				tic=null;
				
				selectedIon=null;
				
				//reset skipping flag and increment counter
				isSkipSpectrumWithBadEncoding = false;
				numSkippedWithBadEncoding++;
				
				return;
			}
			

			if (massArray==null||intensityArray==null) {
				// sometimes SCIEX files skip the binary data and don't populate massArray/intensityArray
				massArray=new double[0];
				intensityArray=new float[0];
			}
			
			if (spectrumRef==null&&msLevel<=1) {
				if (parameters!=null&&parameters.getPrecursorOffsetPPM()!=0.0) {
					double[] deltaArray=General.multiply(massArray, parameters.getPrecursorOffsetPPM()/1000000.0);
					massArray=General.subtract(massArray, deltaArray);
				}
				precursors.add(new PrecursorScan(spectrumName, spectrumIndex, scanStartTime, fraction, scanWindowLowerLimit, scanWindowUpperLimit, ionInjectTime, checkArray(massArray), checkArray(intensityArray), tic));

			} else {
				if (spectrumRef==null) spectrumRef="Unknown";
				
				if (parameters!=null&&parameters.getFragmentOffsetPPM()!=0.0) {
					double[] deltaArray=General.multiply(massArray, parameters.getFragmentOffsetPPM()/1000000.0);
					massArray=General.subtract(massArray, deltaArray);
				}
				
				//Issue 3
				if (isolationWindowTarget==null||isolationWindowLowerOffset==null||isolationWindowUpperOffset==null) {
					if (parameters!=null&&parameters.getPrecursorWindowSize()>0f&&selectedIon!=null) {
						isolationWindowTarget=selectedIon;
						isolationWindowLowerOffset=parameters.getPrecursorWindowSize()/2.0;
						isolationWindowUpperOffset=parameters.getPrecursorWindowSize()/2.0;
					} else {
						Logger.errorLine("Isolation window information missing without precursor window size supplied!");
					}
				}
				
				byte charge=0;
				if (selectedCharge!=null) {
					charge=selectedCharge;
				}
				
				if (parameters!=null&&parameters.getMinIntensity()>0.0f) {
					ArrayList<Peak> peaks=new ArrayList<>();
					for (int i=0; i<intensityArray.length; i++) {
						if (intensityArray[i]>parameters.getMinIntensity()) {
							peaks.add(new Peak(massArray[i], intensityArray[i]));
						}
					}
					Pair<double[], float[]> peakArrays=Peak.toArrays(peaks);
					massArray=peakArrays.x;
					intensityArray=peakArrays.y;
				}
				
				double precursorIsolationMargin = parameters==null?0.0:parameters.getPrecursorIsolationMargin();
				try {
					FragmentScan stripe=new FragmentScan(spectrumName, spectrumRef, spectrumIndex, scanStartTime, fraction, ionInjectTime, isolationWindowTarget-isolationWindowLowerOffset+precursorIsolationMargin, isolationWindowTarget+isolationWindowUpperOffset-precursorIsolationMargin,
							checkArray(massArray), checkArray(intensityArray), charge);
					stripes.add(stripe);
					
					Range range=stripe.getRange();
					TFloatArrayList stripeRTs=retentionTimesByStripe.get(range);
					TFloatArrayList stripeIITs=ionInjectionTimesByStripe.get(range);
					if (stripeRTs==null) {
						stripeRTs=new TFloatArrayList();
						retentionTimesByStripe.put(range, stripeRTs);
						stripeIITs=new TFloatArrayList();
						ionInjectionTimesByStripe.put(range, stripeIITs);
					}
					stripeRTs.add(scanStartTime);
					if (ionInjectTime!=null) {
						stripeIITs.add(ionInjectTime);
					}
				} catch (NullPointerException npe) {
					Logger.errorLine("Potential nulls:");
					Logger.errorLine("spectrumName="+spectrumName);
					Logger.errorLine("spectrumRef="+spectrumRef);
					Logger.errorLine("spectrumIndex="+spectrumIndex);
					Logger.errorLine("scanStartTime="+scanStartTime);
					Logger.errorLine("ionInjectTime="+ionInjectTime);
					Logger.errorLine("isolationWindowTarget="+isolationWindowTarget);
					Logger.errorLine("parameters.getPrecursorIsolationMargin()="+precursorIsolationMargin);
					Logger.errorLine("isolationWindowLowerOffset="+isolationWindowLowerOffset);
					Logger.errorLine("isolationWindowUpperOffset="+isolationWindowUpperOffset);
					Logger.errorLine("massArray="+massArray);
					Logger.errorLine("intensityArray="+intensityArray);
					Logger.errorLine("charge="+charge);
					throw new EncyclopediaException("Null pointer found when creating stripe object", npe);
				}
				
			}

			if (precursors.size()>100||stripes.size()>1000) {
				putBlock(new MSMSBlock(precursors, stripes));
				precursors.clear();
				stripes.clear();
			}

			spectrumName=null;
			spectrumIndex=null;
			spectrumRef=null;

			scanStartTime=null;
			isolationWindowTarget=null;
			isolationWindowLowerOffset=null;
			isolationWindowUpperOffset=null;
			scanWindowLowerLimit=null;
			scanWindowUpperLimit=null;

			compress=false;
			encoding=null;
			isMZ=null;
			dataSB.setLength(0);
			massArray=null;
			intensityArray=null;
			tic=null;
			
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

		if (InstrumentComponent.Type.getTypeByName(qName).isPresent()) {
			instrumentComponentsBuilder.add(currentInstrumentComponentBuilder.build());
			currentInstrumentComponentBuilder = InstrumentComponent.builder();
		} else if ("instrumentConfiguration".equalsIgnoreCase(qName)) {
			instrumentIdToInstrumentComponentBuilder.putAll(currentInstrumentConfigurationIdBuilder.build(), instrumentComponentsBuilder.build());
			currentInstrumentConfigurationIdBuilder = InstrumentId.builder();
			instrumentComponentsBuilder = ImmutableList.builder();
		}

		tagList.remove(tagList.size()-1);
	}
	
	/**
	 * edits in place but still returns the same array for pipelining
	 * @param array
	 * @return
	 */
	private float[] checkArray(float[] array) {
		for (int i = 0; i < array.length; i++) {
			if (!Float.isFinite(array[i])) {
				array[i]=0.0f;
			}
		}
		return array;
	}
	
	/**
	 * edits in place but still returns the same array for pipelining
	 * @param array
	 * @return
	 */
	private double[] checkArray(double[] array) {
		for (int i = 0; i < array.length; i++) {
			if (!Double.isFinite(array[i])) {
				array[i]=0.0f;
			}
		}
		return array;
	}

	@Override
	public void characters(char[] ch, int start, int length) throws SAXException {
		if ("binary".equals(tagList.get(tagList.size()-1))) {
			dataSB.append(ch, start, length);
		}
	}

	@Override
	public void endDocument() throws SAXException {
		putBlock(new MSMSBlock(precursors, stripes));
		putBlock(MSMSBlock.POISON_BLOCK);
		if (numSkippedWithBadEncoding > MAX_BAD_ENCODING_SKIPPED) {
			throw new EncyclopediaException("Skipped too many spectra as a result of bad binary data encoding! File is invalid.");
		} else if (numSkippedWithBadEncoding > 0){
			Logger.errorLine("Skipped " + numSkippedWithBadEncoding + " spectra because of bad binary data encoding.");
		}
	}
	
	@Override
	public void putBlock(MSMSBlock block) {
		try {
			mzmlBlockQueue.put(block);
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
