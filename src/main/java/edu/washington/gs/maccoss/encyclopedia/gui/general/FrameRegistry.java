package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.util.Vector;

import javax.swing.JFrame;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class FrameRegistry {
	BackgroundKeyboardListener keyListener=new BackgroundKeyboardListener();
	Vector<JFrame> frames=new Vector<JFrame>();
	
	public synchronized void registerFrame(JFrame frame) {
		frames.add(frame);
		keyListener.addKeyAndContainerListenerRecursively(frame);
	}
	
	public synchronized void removeFrame(JFrame frame) {
		frames.remove(frame);
		keyListener.removeKeyAndContainerListenerRecursively(frame);
	}
	
	public synchronized void removeFrameAndQuit(JFrame frame) {
		removeFrame(frame);
		
		if (frames.size()==0) {
			Logger.logLine("Last frame closed so quitting.");
			System.exit(0);
		} else {
			Logger.logLine("Frame closed, "+frames.size()+" frames still open.");
		}
	}
	
	public synchronized void removeFrameAndForceQuit(JFrame frame) {
		removeFrame(frame);

		Logger.logLine("Force quitting.");
		System.exit(0);
	}
}
