package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SpectrumPeakFilter;

public class MzmlToDIAConsumer implements Runnable {
	private final BlockingQueue<MzmlBlock> mzmlBlockQueue;
	private final StripeFile stripeFile;
	private final SearchParameters parameters;

	public MzmlToDIAConsumer(BlockingQueue<MzmlBlock> mzmlBlockQueue, StripeFile stripeFile, SearchParameters parameters) {
		this.mzmlBlockQueue=mzmlBlockQueue;
		this.stripeFile=stripeFile;
		this.parameters=parameters;
	}

	@Override
	public void run() {
		try {
			float totalPrecursorTIC=0.0f;
			while (true) {
				MzmlBlock block=mzmlBlockQueue.take();
				if (MzmlBlock.POISON_BLOCK==block) break;
				
				for (PrecursorScan precursor : block.getPrecursors()) {
					totalPrecursorTIC+=precursor.getTIC();
				}
				
				stripeFile.addPrecursor(block.getPrecursors());
				
				ArrayList<Stripe> stripes=block.getStripes();
				if (parameters.isFilterPeaklists()) {
					ArrayList<Stripe> filtered=new ArrayList<>();
					for (Stripe stripe : stripes) {
						filtered.add(SpectrumPeakFilter.filterPeaks(stripe));
					}
					stripes=filtered;
				}
				stripeFile.addStripe(stripes);
			}
			
			stripeFile.addMetadata(StripeFile.TOTAL_PRECURSOR_TIC_ATTRIBUTE, Float.toString(totalPrecursorTIC));
			
		} catch (InterruptedException ie) {
			Logger.errorLine("DIA writing interrupted!");
			Logger.errorException(ie);
		} catch (IOException ioe) {
			throw new EncyclopediaException("DIA writing IO error!", ioe);
		} catch (SQLException sqle) {
			throw new EncyclopediaException("DIA writing SQL error!", sqle);
		}
	}
}
