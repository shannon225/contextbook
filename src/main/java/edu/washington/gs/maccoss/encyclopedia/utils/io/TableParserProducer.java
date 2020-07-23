package edu.washington.gs.maccoss.encyclopedia.utils.io;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

import java.io.*;
import java.util.*;
import java.util.concurrent.BlockingQueue;

public class TableParserProducer implements Runnable {
	public static final HashMap<String, String> POISON_BLOCK=new HashMap<String, String>();
	
	private final BlockingQueue<Map<String, String>> blockingQueue;
	private final BufferedReader in;
	private final String delim;
	private final int numConsumers;
	private final Optional<File> optionalFile;
	
	/**
	 * Note: does not respect either escaping or quotations! 
	 * @param f
	 * @param delim
	 * @return
	 */
	public TableParserProducer(BlockingQueue<Map<String, String>> blockingQueue, File f, String delim, int numConsumers) {
		this.blockingQueue=blockingQueue;
		try {
			this.optionalFile=Optional.of(f);
			in=new BufferedReader(new FileReader(f));
		} catch (IOException ioe) {
			throw new EncyclopediaException(ioe);
		}
		this.delim=delim;
		this.numConsumers=numConsumers;
	}

	/**
	 * mainly for testing purposes!
	 * @param blockingQueue
	 * @param s
	 * @param delim
	 * @param numConsumers
	 */
	public TableParserProducer(BlockingQueue<Map<String, String>> blockingQueue, InputStream s, String delim, int numConsumers) {
		this.blockingQueue=blockingQueue;
		this.in=new BufferedReader(new InputStreamReader(s));
		this.delim=delim;
		this.numConsumers=numConsumers;
		this.optionalFile=Optional.empty();
	}

	public BlockingQueue<Map<String, String>> getQueue() {
		return blockingQueue;
	}

	public void run() {
		
		try {
			// header
			String eachline = in.readLine();

			if (eachline!=null) {
				// trim whitespace from the header before splitting out column names
				eachline = eachline.trim();

				List<String> headers=Arrays.asList(eachline.split(delim, -1));

				while ((eachline=in.readLine())!=null) {
					if (eachline.trim().length()==0) {
						continue;
					}

					HashMap<String, String> row=new HashMap<String, String>();

					List<String> data=Arrays.asList(eachline.split(delim, -1));
					int count=0;
					for (String entry : data) {
						if (count>=headers.size()) {
							break;
						}
						row.put(headers.get(count), entry);
						count++;
					}

					blockingQueue.add(row);
				}
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
