package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SpectrumPeakFilter;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

public class MSMSToDIAConsumer implements Runnable {
	private final BlockingQueue<MSMSBlock> mzmlBlockQueue;
	private final StripeFile stripeFile;
	private final SearchParameters parameters;

	private Throwable error;

	public MSMSToDIAConsumer(BlockingQueue<MSMSBlock> mzmlBlockQueue, StripeFile stripeFile, SearchParameters parameters) {
		this.mzmlBlockQueue=mzmlBlockQueue;
		this.stripeFile=stripeFile;
		this.parameters=parameters;
	}

	@Override
	public void run() {
		try {
			float totalPrecursorTIC=0.0f;
			while (true) {
				MSMSBlock block=mzmlBlockQueue.take();
				if (MSMSBlock.POISON_BLOCK==block) break;

				for (PrecursorScan precursor : block.getPrecursors()) {
					totalPrecursorTIC+=precursor.getTIC();
				}

				stripeFile.addPrecursor(block.getPrecursors());

				ArrayList<FragmentScan> stripes=block.getStripes();
				if (parameters.isFilterPeaklists()) {
					ArrayList<FragmentScan> filtered=new ArrayList<>();
					for (FragmentScan stripe : stripes) {
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
		} catch (Throwable t) {
			Logger.errorLine("Error writing to DIA!");
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
}
