package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.SearchTestSupport;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursorWithProteins;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.dia.FragmentIonConsistencyCharter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import junit.framework.TestCase;

public class PeptidePropertiesModelTest extends TestCase {
	public static void main(String[] args) throws Exception {
		PeptidePropertiesModel model=PeptidePropertiesModel.getModel();
		
		File libraryFile=new File("/Users/searleb/Downloads/pan_human_library.dlib");
		LibraryFile library=new LibraryFile();
		library.openFile(libraryFile);

    	SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
    	AminoAcidConstants aaConstants=parameters.getAAConstants();
    	ArrayList<LibraryEntry> entries=library.getAllEntries(false, aaConstants);
    	
    	ArrayList<PeptidePrecursorWithProteins> peptidesWithProteins=new ArrayList<PeptidePrecursorWithProteins>();
        for (LibraryEntry peptide : entries) {
        	if (peptide.getPeptideSeq().length()>PeptideEncoding.MAX_PEPTIDE_LENGTH-2||peptide.getPrecursorCharge()>PeptideEncoding.MAX_CHARGE) {
        		continue;
        	}
        	peptidesWithProteins.add(peptide);
		}
    	ArrayList<LibraryEntry> predicted=model.predict(peptidesWithProteins, aaConstants);

    	int index = 30;
		LibraryEntry top=new AnnotatedLibraryEntry((LibraryEntry)peptidesWithProteins.get(index), parameters);
    	LibraryEntry bottom=new AnnotatedLibraryEntry(predicted.get(index), parameters);
    	System.out.println(bottom.getPeptideModSeq()+", +"+bottom.getPrecursorCharge());
    	Charter.launchChart(new AnnotatedLibraryEntry(FragmentIonConsistencyCharter.getButterfly(top, bottom), parameters, true));
    	
    	System.out.println(peptidesWithProteins.size()+"/"+predicted.size());
    	
    	ArrayList<XYPoint> rts=new ArrayList<XYPoint>();
    	for (int i = 0; i < peptidesWithProteins.size(); i++) {
			rts.add(new XYPoint(((LibraryEntry)peptidesWithProteins.get(i)).getRetentionTime()/60f, predicted.get(i).getRetentionTime()/60f));
		}

		RetentionTimeFilter filter=RetentionTimeFilter.getFilter(rts);
		filter.plot(rts, Optional.ofNullable((File)null));
	}

	public void testSmoke() throws Exception {
		PeptidePropertiesModel model=PeptidePropertiesModel.getModel();
		LibraryFile library=SearchTestSupport.getResultLibrary();

    	SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
    	AminoAcidConstants aaConstants=parameters.getAAConstants();
    	ArrayList<LibraryEntry> entries=library.getAllEntries(false, aaConstants);
    	
    	ArrayList<PeptidePrecursorWithProteins> peptidesWithProteins=new ArrayList<PeptidePrecursorWithProteins>(entries);
    	ArrayList<LibraryEntry> predicted=model.predict(peptidesWithProteins, aaConstants);

    	assert(entries.size()==predicted.size());
	}
}
