package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlBlock;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class OverlapDeconvoluter implements Runnable {
	private final BlockingQueue<MzmlBlock> inputQueue;
	private final BlockingQueue<MzmlBlock> outputQueue;

	public OverlapDeconvoluter(BlockingQueue<MzmlBlock> inputQueue, BlockingQueue<MzmlBlock> outputQueue) {
		this.inputQueue=inputQueue;
		this.outputQueue=outputQueue;
	}

	public void run() {
		try {
			while (true) {
				MzmlBlock block=inputQueue.take();
				if (MzmlBlock.POISON_BLOCK==block) break;
				// FIXME do work
				outputQueue.put(block);
			}
		} catch (InterruptedException ie) {
			Logger.errorLine("DIA writing interrupted!");
			Logger.errorException(ie);
		}
	}
	
}
