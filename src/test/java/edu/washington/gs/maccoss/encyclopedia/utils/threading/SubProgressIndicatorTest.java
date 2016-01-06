package edu.washington.gs.maccoss.encyclopedia.utils.threading;

import junit.framework.TestCase;

public class SubProgressIndicatorTest extends TestCase {
	public void testSubProgress() {
		EmptyProgressIndicator progress=new EmptyProgressIndicator();
		SubProgressIndicator sp1=new SubProgressIndicator(progress, 0.2f);
		sp1.update("", 0.0f); assertEquals(progress.getTotalProgress(),sp1.getTotalProgress());assertEquals(0.0f,sp1.getTotalProgress(), 0.00001f);
		sp1.update("", 0.5f); assertEquals(progress.getTotalProgress(),sp1.getTotalProgress());assertEquals(0.1f,sp1.getTotalProgress(), 0.00001f);
		sp1.update("", 0.9f); assertEquals(progress.getTotalProgress(),sp1.getTotalProgress());assertEquals(0.18f,sp1.getTotalProgress(), 0.00001f);
		sp1.update("", 1.0f); assertEquals(progress.getTotalProgress(),sp1.getTotalProgress());assertEquals(0.2f,sp1.getTotalProgress(), 0.00001f);
		

		SubProgressIndicator sp2=new SubProgressIndicator(progress, 0.2f);
		sp2.update("", 0.0f); assertEquals(progress.getTotalProgress(),sp2.getTotalProgress());assertEquals(0.2f,sp2.getTotalProgress(), 0.00001f);
		sp2.update("", 0.5f); assertEquals(progress.getTotalProgress(),sp2.getTotalProgress());assertEquals(0.3f,sp2.getTotalProgress(), 0.00001f);
		sp2.update("", 0.9f); assertEquals(progress.getTotalProgress(),sp2.getTotalProgress());assertEquals(0.38f,sp2.getTotalProgress(), 0.00001f);
		sp2.update("", 1.0f); assertEquals(progress.getTotalProgress(),sp2.getTotalProgress());assertEquals(0.4f,sp2.getTotalProgress(), 0.00001f);

		SubProgressIndicator sp3=new SubProgressIndicator(progress, 0.2f);
		sp3.update("", 0.0f); assertEquals(progress.getTotalProgress(),sp3.getTotalProgress());assertEquals(0.4f,sp3.getTotalProgress(), 0.00001f);
		sp3.update("", 0.5f); assertEquals(progress.getTotalProgress(),sp3.getTotalProgress());assertEquals(0.5f,sp3.getTotalProgress(), 0.00001f);
		
		SubProgressIndicator sp4=new SubProgressIndicator(sp3, 0.5f);
		sp4.update("", 0.5f); assertEquals(progress.getTotalProgress(),sp4.getTotalProgress());assertEquals(0.55f,sp4.getTotalProgress(), 0.00001f);
		sp4.update("", 1.0f); assertEquals(progress.getTotalProgress(),sp4.getTotalProgress());assertEquals(0.6f,sp4.getTotalProgress(), 0.00001f);
	}

}
