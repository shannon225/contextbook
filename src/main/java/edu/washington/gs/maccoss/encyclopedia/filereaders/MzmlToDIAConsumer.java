package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.IOException;
import java.sql.SQLException;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class MzmlToDIAConsumer implements Runnable {
	private final BlockingQueue<MzmlBlock> mzmlBlockQueue;
	private final StripeFile stripeFile;

	public MzmlToDIAConsumer(BlockingQueue<MzmlBlock> mzmlBlockQueue, StripeFile stripeFile) {
		this.mzmlBlockQueue=mzmlBlockQueue;
		this.stripeFile=stripeFile;
	}

	@Override
	public void run() {
		try {
			while (true) {
				MzmlBlock block=mzmlBlockQueue.take();
				if (MzmlBlock.POISON_BLOCK==block) break;
				
				stripeFile.addPrecursor(block.getPrecursors());
				stripeFile.addStripe(block.getStripes());
			}
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
