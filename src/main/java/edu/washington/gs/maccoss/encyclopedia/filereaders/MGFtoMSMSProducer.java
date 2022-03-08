package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.StringTokenizer;
import java.util.concurrent.BlockingQueue;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.io.ProgressInputStream;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;

public class MGFtoMSMSProducer implements MSMSProducer {

	private final File mgfFile;
	private final BlockingQueue<MSMSBlock> mgfBlockQueue;
	private Throwable error;
	
	public MGFtoMSMSProducer(File mgfFile, BlockingQueue<MSMSBlock> mzmlBlockQueue) {
		this.mgfFile=mgfFile;
		this.mgfBlockQueue=mzmlBlockQueue;
	}

	@Override
	public void run() {
		try {
			final ProgressInputStream stream = new ProgressInputStream(new FileInputStream(mgfFile));
			final long length = mgfFile.length();

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

			BufferedReader in=new BufferedReader(new InputStreamReader(stream));
			String eachline;
			String title=null;
			float rtInSecs=0.0f;
			double mz=0.0;
			byte charge=0;
			ArrayList<Peak> peaks=new ArrayList<Peak>();
			ArrayList<FragmentScan> spectra=new ArrayList<FragmentScan>();
			
			int count=0;
			while ((eachline=in.readLine())!=null) {
				if (eachline.equals("BEGIN IONS")) {
					if (peaks.size()>0) {
						// add spectrum
						count++;
						Pair<double[], float[]> arrays=Peak.toArrays(peaks);
						FragmentScan scan=new FragmentScan(title, "", count, rtInSecs, 1, null, mz-0.01, mz+0.01, arrays.x, arrays.y, charge);
						spectra.add(scan);
						
						if (spectra.size()>1000) {
							putBlock(new MSMSBlock(new ArrayList<>(), spectra));
							spectra.clear();
						}
					} else if (title!=null) {
						Logger.errorLine("No peaks assigned to ["+title+']');
					}
					peaks.clear();
				}
				if (eachline.startsWith("TITLE=")) {
					title=eachline.substring("TITLE=".length());
				} else if (eachline.startsWith("RTINSECONDS=")) {
					rtInSecs=Float.parseFloat(eachline.substring("RTINSECONDS=".length()));
					
				} else if (eachline.startsWith("PEPMASS=")) {
					mz=Double.parseDouble(eachline.substring("RTINSECONDS=".length(), eachline.indexOf(' ')));
					
				} else if (eachline.startsWith("CHARGE=")) {
					if (eachline.endsWith("+")) {
					charge=Byte.parseByte(eachline.substring("CHARGE=".length(), eachline.indexOf('+')));
					} else if (eachline.endsWith("-")) {
						charge=Byte.parseByte(eachline.substring("CHARGE=".length(), eachline.indexOf('-')));
					} else {
						charge=Byte.parseByte(eachline.substring("CHARGE=".length()));
					}
					
				} else if (Character.isDigit(eachline.charAt(0))) {
					StringTokenizer st=new StringTokenizer(eachline);
					double mass=Double.parseDouble(st.nextToken());
					float intensity=Float.parseFloat(st.nextToken());
					peaks.add(new Peak(mass, intensity));
				}
			}
			if (spectra.size()>0) {
				putBlock(new MSMSBlock(new ArrayList<>(), spectra));
			}
			
			stream.close();

			putBlock(MSMSBlock.POISON_BLOCK);
		} catch (Throwable t) {
			Logger.errorLine("MGF reading failed!");
			Logger.errorException(t);

			this.error = t;
		}
	}
	
	@Override
	public void putBlock(MSMSBlock block) {
		try {
			mgfBlockQueue.put(block);
		} catch (InterruptedException ie) {
			Logger.errorLine("MGF reading interrupted!");
			Logger.errorException(ie);
		}
	}

	public boolean hadError() {
		return null != error;
	}

	public Throwable getError() {
		return error;
	}

}
