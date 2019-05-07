package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.awt.Component;
import java.text.DecimalFormat;

import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.table.DefaultTableCellRenderer;

public class TICTableCellRenderer extends DefaultTableCellRenderer {
	private static final long serialVersionUID = 1L;
	private static final DecimalFormat decFormat = new DecimalFormat("#.00");

	public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus,
			int row, int column) {
		super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
		this.setHorizontalAlignment(SwingConstants.RIGHT);
		
		if (value instanceof Number) {
			String text = toScientific(((Number)value).doubleValue());
			setText(text);
			setToolTipText(text);
		} else if (value != null) {
			setText(value.toString());
			setToolTipText((String) value);
		} else {
			setText("???");
			setToolTipText("null value");
		}

		return this;
	}

	private static String toScientific(double num) {
		if (num==0.0) return "0.00e0";
		String result = "";

		int power = (int)(Math.log(num)/Math.log(10));
		if (power < 0) power--;
		double fraction = num/Math.pow(10, power);
		result += decFormat.format(fraction)+"e"+power;
		return result;
	}
}