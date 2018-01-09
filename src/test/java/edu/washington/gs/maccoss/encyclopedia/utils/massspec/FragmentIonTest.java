package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import junit.framework.TestCase;

public class FragmentIonTest extends TestCase {
	public void testArchive() {
		FragmentIon[] ions=new FragmentIon[] {
			new FragmentIon(100.0, (byte)1, IonType.b),
			new FragmentIon(200.0, (byte)2, IonType.y),
			new FragmentIon(300.0, (byte)3, IonType.bp2),
			new FragmentIon(400.0, (byte)4, IonType.yNL),
			new FragmentIon(500.0, (byte)5, IonType.bp2NL),	
		};
		
		String s=FragmentIon.toArchiveString(ions);
		
		FragmentIon[] extracted=FragmentIon.fromArchiveString(s);
		assertEquals(ions.length, extracted.length);
		for (int i=0; i<extracted.length; i++) {
			assertEquals(ions[i], extracted[i]);
		}
	}
	
	public void testUniqueModificationIons() {
		String peptideModSeq="EM[+15.994915]DEAATAEER";
		byte charge=2;
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		
		FragmentationModel model=PeptideUtils.getPeptideModel(peptideModSeq, parameters.getAAConstants());
		FragmentIon[] ions=model.getPrimaryIonObjects(parameters.getFragType(), charge, true);
		assertEquals(22, ions.length);

		Range precursorIsolationRange=new Range(628.535583496094f, 640.541076660156f);
		Optional<FragmentIon[]> modificationSpecificIons=model.getModificationSpecificIonObjects(precursorIsolationRange, parameters.getFragType(), charge, true);
		assertFalse(modificationSpecificIons.isPresent());

		precursorIsolationRange=new Range(616.535583496094f, 640.541076660156f);
		modificationSpecificIons=model.getModificationSpecificIonObjects(precursorIsolationRange, parameters.getFragType(), charge, true);
		assertTrue(modificationSpecificIons.isPresent());
		assertEquals(12, modificationSpecificIons.get().length);
	}
}
