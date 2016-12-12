package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

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
		for (int i=0; i<extracted.length; i++) {
			assertEquals(ions[i], extracted[i]);
		}
	}
}
