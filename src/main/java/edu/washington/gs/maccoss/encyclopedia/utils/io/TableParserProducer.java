package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class TableParserProducer implements Runnable {
	public static final HashMap<String, String> POISON_BLOCK=new HashMap<String, String>();
	
	private final BlockingQueue<Map<String, String>> blockingQueue;
	private final File f;
	private final String delim;
	private final int numConsumers;
	
	/**
	 * Note: does not respect either escaping or quotations! 
	 * @param f
	 * @param delim
	 * @return
	 */
	public TableParserProducer(BlockingQueue<Map<String, String>> blockingQueue, File f, String delim, int numConsumers) {
		this.blockingQueue=blockingQueue;
		this.f=f;
		this.delim=delim;
		this.numConsumers=numConsumers;
	}
	
	public void run() {
		BufferedReader in=null;
		
		try {
			in=new BufferedReader(new FileReader(f));
			
			String eachline=in.readLine(); // header
			
			if (eachline!=null) {
				StringTokenizer st=new StringTokenizer(eachline, delim);
				ArrayList<String> headers=new ArrayList<String>();
				while (st.hasMoreTokens()) {
					headers.add(st.nextToken());
				}

				while ((eachline=in.readLine())!=null) {
					if (eachline.trim().length()==0) {
						continue;
					}

					HashMap<String, String> row=new HashMap<String, String>();

					st=new StringTokenizer(eachline, delim);
					int count=0;
					while (st.hasMoreTokens()) {
						String entry=st.nextToken();
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
			Logger.errorLine("I/O Error found reading table ["+f.getName()+"]");
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
