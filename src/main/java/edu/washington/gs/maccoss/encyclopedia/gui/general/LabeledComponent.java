package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.BorderLayout;
import java.awt.Color;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;

public class LabeledComponent extends JPanel {
	private static final long serialVersionUID=1L;

	public LabeledComponent(String label, JComponent c) {
		super(new BorderLayout());

		add(new JLabel("<html><p style=\"font-size:10px; font-family: Helvetica, sans-serif\">"+label+":"), BorderLayout.WEST);
		add(c, BorderLayout.CENTER);
		this.setOpaque(true);
		this.setBackground(Color.WHITE);
	}

}
