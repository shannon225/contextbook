package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class FileConcatenator {
	private static final String concatenationBoundary="~~~~FILEBREAK~~~~";
	
	public static int getNumberOfSubFiles(File concatFile) throws IOException {
		BufferedReader in=new BufferedReader(new FileReader(concatFile));
		int count=1;
		String eachline;
		while ((eachline=in.readLine())!=null) {
			if (eachline.equals(concatenationBoundary)) {
				count++;
			}
		}
		
		in.close();
		return count;
	}

	public static File extractFile(File concatFile, int n) throws IOException {
		File tempFile=File.createTempFile("encyclopedia_single_", ".single");
		tempFile.deleteOnExit();
		extractFile(concatFile, n, tempFile);
		return tempFile;
	}
	
	public static void extractFile(File concatFile, int n, File singleFile) throws IOException {
		BufferedWriter out=new BufferedWriter(new FileWriter(singleFile));
		BufferedReader in=new BufferedReader(new FileReader(concatFile));
		int count=1;
		String eachline;
		while ((eachline=in.readLine())!=null) {
			if (eachline.equals(concatenationBoundary)) {
				count++;
				if (count>n) {
					break;
				} else {
					continue;
				}
			}
			
			if (count==n) {
				out.write(eachline);
				out.newLine();
			}
		}
		
		in.close();
		out.flush();
		out.close();
	}

	public static void saveConcatenatedFile(File concatFile, File...singleFiles) throws IOException {
		int count=0;
		if (concatFile.exists()&&concatFile.length()>0) {
			// append to file
			count=1;
		}
		BufferedWriter out=new BufferedWriter(new FileWriter(concatFile, true));
		for (File file : singleFiles) {
			if (count>0) {
				out.write(concatenationBoundary);
				out.newLine();
			}
			
			BufferedReader in=new BufferedReader(new FileReader(file));
			String eachline;
			while ((eachline=in.readLine())!=null) {
				out.write(eachline);
				out.newLine();
			}
			in.close();
			
			count++;
		}
		out.flush();
		out.close();
	}
}
