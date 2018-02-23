package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.*;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class LibraryBackground implements LibraryBackgroundInterface {
	private final int[] background=new int[4000];
	private final int total;
	
	public static void main(String[] args) throws Exception {
		File libraryFile=new File("/Users/searleb/Documents/school/encyclopedia_manuscript/22oct2017_hela_serum_timecourse_narrow_library.elib");
		LibraryFile file=new LibraryFile();
		file.openFile(libraryFile);
		ArrayList<LibraryEntry> entries=file.getEntries(new Range(700.57f, 712.57f), false, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()));
		LibraryBackground background=new LibraryBackground(entries);
		for (int i=100; i<=1000; i++) {
			System.out.println(i+"\t"+(1.0f/background.getFraction(i+0.1))); // to avoid rounding errors
		}
	}

	public LibraryBackground(Collection<String> peptides, PecanSearchParameters params) {
		Arrays.fill(background, 1); // add a pseudocount
		for (String sequence : peptides) {
			for (byte charge=params.getMinCharge(); charge<=params.getMaxCharge(); charge++) {
				FragmentationModel model=PeptideUtils.getPeptideModel(sequence, params.getAAConstants());
				double[] ions=model.getPrimaryIons(params.getFragType(), charge, false);
				for (double ion : ions) {
					int index=(int)ion; // truncate
					if (index<background.length) {
						background[index]++;
					}
				}
			}
		}
		int t=0;
		for (int i=0; i<background.length; i++) {
			t+=background[i];
		}
		this.total=t;
	}

	public LibraryBackground(ArrayList<LibraryEntry> entries) {
		Arrays.fill(background, 1); // add a pseudocount
		for (LibraryEntry entry : entries) {
			double[] masses=entry.getMassArray();
			for (double mass : masses) {
				int index=(int)mass; // truncate
				if (index<background.length) {
					background[index]++;
				}
			}
		}
		int t=0;
		for (int i=0; i<background.length; i++) {
			t+=background[i];
		}
		this.total=t;
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackgroundInterface#getFraction(double)
	 */
	@Override
	public float getFraction(double mass) {
		int index=(int)mass; // truncate
		int count=index>=background.length?1:background[index];
		return (total/(float)background.length)/(float)count;
	}
}
