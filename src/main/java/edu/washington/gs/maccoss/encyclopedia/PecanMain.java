package edu.washington.gs.maccoss.encyclopedia;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.JFrame;

import edu.washington.gs.maccoss.encyclopedia.gui.pecan.PecanPanel;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class PecanMain {
	public static void main(String[] args) {
		final JFrame f=new JFrame("Pecan Peptide Centric Analysis");
		f.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});

		f.getContentPane().add(new PecanPanel(), BorderLayout.CENTER);
		
		f.pack();
		f.setSize(new Dimension(1100, 800));
		f.setVisible(true);

		Logger.logLine("Pecanpie Graphical Interface");
	}

}