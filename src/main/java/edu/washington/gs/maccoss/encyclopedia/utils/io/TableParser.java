package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class TableParser {
	/**
	 * Note: does not respect either escaping or quotations! 
	 * @param f
	 * @param delim
	 * @return
	 */
	public static void readTable(File f, String delim, TableParserMuscle muscle) {
		BufferedReader in=null;
		
		try {
			in=new BufferedReader(new FileReader(f));
			
			String eachline=in.readLine(); // header
			StringTokenizer st=new StringTokenizer(eachline, delim);
			ArrayList<String> headers=new ArrayList<String>();
			while (st.hasMoreTokens()) {
				headers.add(st.nextToken());
			}
			
			while ((eachline=in.readLine())!=null) {
				if (eachline.trim().length()==0) {
					continue;
				}

				HashMap<String, String> map=new HashMap<String, String>();
				
				st=new StringTokenizer(eachline, delim);
				int count=0;
				while (st.hasMoreTokens()) {
					String entry=st.nextToken();
					if (count>=headers.size()) {
						break;
					}
					map.put(headers.get(count), entry);
					count++;
					
				}
				
				muscle.processRow(map);
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
