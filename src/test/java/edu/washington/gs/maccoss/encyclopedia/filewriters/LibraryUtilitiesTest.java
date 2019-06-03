package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;

public class LibraryUtilitiesTest {
	public static void main2(String[] args) throws Exception {

		File twoDLC=new File("/Users/searleb/Downloads/msms.dlib");

		LibraryFile library=new LibraryFile();
		library.openFile(twoDLC);

		float[] timeBoundaries=new float[] {0f, 2f, 32f, 62f, 92f, 122f, 152f, 182f, 212f, 242f, 272f, 302f, 332f, 362f, 392f, 1000f};
		for (int i=1; i<timeBoundaries.length; i++) {
			File saveFile=new File(twoDLC.getParentFile(), "msms_"+Math.round(timeBoundaries[i-1])+"to"+Math.round(timeBoundaries[i])+".dlib");
			System.out.println("writing "+saveFile.getName());
			LibraryUtilities.subsetLibrary(saveFile, timeBoundaries[i-1]*60f, timeBoundaries[i]*60f, new HashSet<>(), library);
		}

		library.close();
	}
}
