package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.util.ArrayList;

import javax.swing.table.AbstractTableModel;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;

public class DIAScanTableModel extends AbstractTableModel {
	private static final long serialVersionUID=1L;
	
	private final String[] columns=new String[] {"#", "Spectrum Name", "Scan Start Time", "Precursor M/Z"};

	ArrayList<Spectrum> entries=new ArrayList<Spectrum>();
	
	public void updateEntries(ArrayList<Spectrum> newEntries) {
		entries.clear();
		entries.addAll(newEntries);
		fireTableDataChanged();
	}
	
	public Spectrum getSelectedRow(int rowIndex) {
		return entries.get(rowIndex);
	}

	@Override
	public int getRowCount() {
		return entries.size();
	}

	@Override
	public int getColumnCount() {
		return columns.length;
	}
	
	@Override
	public String getColumnName(int column) {
		return columns[column];
	}
	
	@Override
	public Class<?> getColumnClass(int columnIndex) {
		switch (columnIndex) {
			case 0: return Integer.class;
			case 1: return String.class;
			case 2: return Float.class;
			case 3: return Double.class;
		}
		return Object.class;
	}

	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		Spectrum entry=getSelectedRow(rowIndex);
		
		switch (columnIndex) {
			case 0: return rowIndex;
			case 1: return entry.getSpectrumName();
			case 2: return entry.getScanStartTime()/60f;
			case 3: return entry.getPrecursorMZ();
		}
		return null;
	}
	
	
}
