package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.BlockingQueue;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class DIAtoMSMSProducer implements MSMSProducer {
	private static final int STRIPE_CHUNKS=20;

	private final File diaFile;
	private StripeFile stripeFile; 
	private final BlockingQueue<MSMSBlock> outputBlockQueue;
	private final HashMap<String, String> metadata=new HashMap<String, String>();

	private Throwable error;

	public DIAtoMSMSProducer(File diaFile, BlockingQueue<MSMSBlock> outputBlockQueue) {
		this.diaFile=diaFile;
		this.outputBlockQueue=outputBlockQueue;
	}

	@Override
	public void run() {
		try {
			stripeFile=new StripeFile();
			stripeFile.openFile(diaFile);

			float totalRTInSec=stripeFile.getGradientLength();
			
			float start=-Float.MAX_VALUE;
			float stop=-Float.MAX_VALUE;
			for (int i = 1; i < STRIPE_CHUNKS; i++) {
				start=Float.intBitsToFloat(Float.floatToIntBits(stop)+1);
				stop=totalRTInSec*(i/(float)STRIPE_CHUNKS);
				processStripeFile(stripeFile, start, stop);
			}
			start=Float.intBitsToFloat(Float.floatToIntBits(stop)+1);
			processStripeFile(stripeFile, start, Float.MAX_VALUE);
			
			metadata.putAll(stripeFile.getMetadata());
			
			stripeFile.close();
			
		} catch (IOException ioe) {
			error=ioe;
			throw new EncyclopediaException("DIA writing IO error!", ioe);
		} catch (SQLException sqle) {
			error=sqle;
			throw new EncyclopediaException("DIA writing SQL error!", sqle);
		} catch (DataFormatException dfe) {
			error=dfe;
			throw new EncyclopediaException("DIA writing DFE error!", dfe);
		}
	}

	public boolean hadError() {
		return null != error;
	}

	public Throwable getError() {
		return error;
	}
	
	public HashMap<String, String> getMetadata() {
		return metadata;
	}
	
	private void processStripeFile(StripeFileInterface stripeFile, float minRT, float maxRT) throws IOException, SQLException, DataFormatException {
		ArrayList<PrecursorScan> precursors=stripeFile.getPrecursors(minRT, maxRT);
		ArrayList<FragmentScan> fragments=stripeFile.getStripes(new Range(-Float.MAX_VALUE, Float.MAX_VALUE), minRT, maxRT, false);

		putBlock(new MSMSBlock(precursors, fragments));
	}

	@Override
	public void putBlock(MSMSBlock block) {
		try {
			outputBlockQueue.put(block);
		} catch (InterruptedException ie) {
			Logger.errorLine("Mzml reading interrupted!");
			Logger.errorException(ie);
		}
	}
}
