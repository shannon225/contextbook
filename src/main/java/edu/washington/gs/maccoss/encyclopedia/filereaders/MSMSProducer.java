package edu.washington.gs.maccoss.encyclopedia.filereaders;

public interface MSMSProducer extends Runnable {
	public void putBlock(MSMSBlock block);
}
