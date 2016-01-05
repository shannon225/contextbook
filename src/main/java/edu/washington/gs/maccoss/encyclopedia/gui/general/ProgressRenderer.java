package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.Color;
import java.awt.Component;

import javax.swing.BorderFactory;
import javax.swing.JProgressBar;
import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.table.DefaultTableCellRenderer;

import edu.washington.gs.maccoss.encyclopedia.utils.io.ProgressMessage;

public class ProgressRenderer extends DefaultTableCellRenderer {
	private static final long serialVersionUID=1L;

	private static final int REFINEMENT=1000;
	
	private final JProgressBar bar=new JProgressBar(0, REFINEMENT);

	public ProgressRenderer() {
		super();
		setOpaque(true);
		setHorizontalAlignment(SwingConstants.CENTER);
		
		bar.setBorder(BorderFactory.createEmptyBorder(1, 1, 1, 1));
		bar.setStringPainted(true);
	}

	@Override
	public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
		ProgressMessage f=(ProgressMessage)value;
		String text=f.getMessage();
		if (f.isError()) {
			this.setBackground(Color.pink);
			this.setOpaque(true);
		} else if (f.isFinished()) {
			this.setBackground(new Color(150, 255, 150));
			this.setOpaque(true);
		} else {
			bar.setValue(Math.round(f.getProgress()*REFINEMENT));
			bar.setString(f.getMessage());
			return bar;
		}
		super.getTableCellRendererComponent(table, text, isSelected, hasFocus, row, column);
		return this;
	}
}
