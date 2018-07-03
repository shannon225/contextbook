package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Optional;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class LineParserProducer implements Runnable {
	public static final String POISON_BLOCK="LINE_PARSER_POISON_BLOCK";
	
	private final BlockingQueue<String> blockingQueue;
	private final BufferedReader in;
	private final int numConsumers;
	private final Optional<File> optionalFile;
	
	/**
	 * Note: does not respect either escaping or quotations! 
	 * @param f
	 * @param delim
	 * @return
	 */
	public LineParserProducer(BlockingQueue<String> blockingQueue, File f, int numConsumers) {
		this.blockingQueue=blockingQueue;
		try {
			this.optionalFile=Optional.of(f);
			in=new BufferedReader(new FileReader(f));
		} catch (IOException ioe) {
			throw new EncyclopediaException(ioe);
		}
		this.numConsumers=numConsumers;
	}

	/**
	 * mainly for testing purposes!
	 * @param blockingQueue
	 * @param s
	 * @param delim
	 * @param numConsumers
	 */
	public LineParserProducer(BlockingQueue<String> blockingQueue, InputStream s, int numConsumers) {
		this.blockingQueue=blockingQueue;
		this.in=new BufferedReader(new InputStreamReader(s));
		this.numConsumers=numConsumers;
		this.optionalFile=Optional.empty();
	}
	
	public void run() {
		
		try {
			String eachline;
			while ((eachline=in.readLine())!=null) {
				if (eachline.trim().length()==0) {
					continue;
				}
				blockingQueue.add(eachline);
			}
			
			for (int i=0; i<numConsumers; i++) {
				blockingQueue.add(POISON_BLOCK);
			}

		} catch (IOException ioe) {
			if (optionalFile.isPresent()) {
				Logger.errorLine("I/O Error found reading table ["+optionalFile.get().getName()+"]");
			}
			Logger.errorException(ioe);
		} finally {
			if (in!=null) {
				try {
					in.close();
				} catch (IOException ioe) {
					ioe.printStackTrace();
				}
			}
		}
	}
}
