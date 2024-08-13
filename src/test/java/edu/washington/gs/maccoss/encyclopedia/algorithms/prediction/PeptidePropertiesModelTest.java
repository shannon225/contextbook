package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Optional;

import org.nd4j.linalg.api.ndarray.INDArray;

import edu.washington.gs.maccoss.encyclopedia.algorithms.SearchTestSupport;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaTwoAlignmentLibraryFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursorWithProteins;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.dia.FragmentIonConsistencyCharter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Correlation;
import edu.washington.gs.maccoss.encyclopedia.utils.math.PivotTableGenerator;
import gnu.trove.list.array.TFloatArrayList;
import junit.framework.TestCase;

public class PeptidePropertiesModelTest extends TestCase {
	public static void main(String[] args) throws Exception {
		PeptidePropertiesModel model=PeptidePropertiesModel.getModel();
		File libraryFile=EncyclopediaTwoAlignmentLibraryFactory.getPreAlignmentFile(DigestionEnzyme.getEnzyme("Trypsin")).get();
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

    	int index = 32;
		LibraryEntry top=new AnnotatedLibraryEntry((LibraryEntry)peptidesWithProteins.get(index), parameters);
    	LibraryEntry bottom=new AnnotatedLibraryEntry(predicted.get(index), parameters);
    	System.out.println(bottom.getPeptideModSeq()+", +"+bottom.getPrecursorCharge());
    	Charter.launchChart(new AnnotatedLibraryEntry(FragmentIonConsistencyCharter.getButterfly(top, bottom), parameters, true));
    	
    	System.out.println(peptidesWithProteins.size()+"/"+predicted.size());
    	
    	ArrayList<XYPoint> rts=new ArrayList<XYPoint>();
    	TFloatArrayList spectralAngles=new TFloatArrayList();
    	for (int i = 0; i < peptidesWithProteins.size(); i++) {
			rts.add(new XYPoint(predicted.get(i).getRetentionTime()/60f, ((LibraryEntry)peptidesWithProteins.get(i)).getRetentionTime()/60f));
			

			AnnotatedLibraryEntry acquired=AnnotatedLibraryEntry.getAnnotationsOnly((LibraryEntry)peptidesWithProteins.get(i), parameters);
			spectralAngles.add((float)Correlation.getSpectralAngle(predicted.get(i), acquired, parameters.getFragmentTolerance()));
		}

		RetentionTimeFilter filter=RetentionTimeFilter.getFilter(rts);
		filter.plot(rts, Optional.ofNullable((File)null));
		
		ArrayList<XYPoint> pivot=PivotTableGenerator.createPivotTable(spectralAngles.toArray(), -1f, 1f, 0.05f);
		Charter.launchChart("Spectral Contrast Angle", "Number Of Peptides", false, new XYTrace(pivot, GraphType.area, "Spectral Contrast"));
	}
	

	
	public void testEncoding() throws Exception {
		PeptidePropertiesModel model=PeptidePropertiesModel.getModel();
		LibraryFile library=SearchTestSupport.getResultLibrary();

    	SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
    	AminoAcidConstants aaConstants=parameters.getAAConstants();
    	ArrayList<LibraryEntry> entries=library.getAllEntries(false, aaConstants);
    	
    	ArrayList<PeptidePrecursorWithProteins> peptidesWithProteins=new ArrayList<PeptidePrecursorWithProteins>(entries);
    	ArrayList<LibraryEntry> predicted=model.predict(peptidesWithProteins, aaConstants);
    	
    	for (int i = 0; i < predicted.size(); i++) {
    		
			PeptideEncoding encoding=new PeptideEncoding(entries.get(i), entries.get(i).getScanStartTime(), parameters);
		}
		for (LibraryEntry entry : entries) {
			PeptideEncoding encoding=new PeptideEncoding(entry, entry.getScanStartTime(), parameters);
			INDArray encodeInput = encoding.encodeInput();
			INDArray encodeResult = encoding.encodeResult();
			
//			System.out.println("Input: ");
//			System.out.println(General.toString(encodeInput.toFloatVector()));
//
//			System.out.println("\nOutput: ");
//			System.out.println(encodeResult);
			
			assertEquals(1.0f, encoding.score(encodeResult));
			
			LibraryEntry result=encoding.outputToEntry(new HashSet<String>(), encodeResult, aaConstants);
			System.out.println(Correlation.getSpectralAngle(entry, result, new MassTolerance(10)));
			
		}
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
