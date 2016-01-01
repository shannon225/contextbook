package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.Font;

import javax.swing.JLabel;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class MemoryMonitor extends JLabel {
	private static final long serialVersionUID=1L;
	
	private static final int mb=1024*1024;
	
	public MemoryMonitor() {
		setHorizontalAlignment(SwingConstants.RIGHT);
		setFont(new Font("sansserif", Font.PLAIN, 10));
		setText(toString());
	}

	public void start() {
		Thread t=new Thread(new Runnable() {
			@Override
			public void run() {
				while (true) {
					try {
						Thread.sleep(500);
					} catch (InterruptedException ie) {
						Logger.errorLine("Error monitoring memory!");
						Logger.errorException(ie);
					}

					SwingUtilities.invokeLater(new Runnable() {
						@Override
						public void run() {
							setText(MemoryMonitor.this.toString());
						}
					});
				}
			}
		});
		t.setDaemon(true);
		t.start();
	}
	
	public String toString() {
		Runtime instance=Runtime.getRuntime();
		long usedMemory=instance.totalMemory()-instance.freeMemory();
		String text=(usedMemory/mb)+" of "+(instance.maxMemory()/mb)+" MB used";
		return text;
	}
}
