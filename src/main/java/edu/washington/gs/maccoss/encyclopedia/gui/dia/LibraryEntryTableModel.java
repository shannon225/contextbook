package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.util.ArrayList;

import javax.swing.table.AbstractTableModel;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;

public class LibraryEntryTableModel extends AbstractTableModel {
	private static final long serialVersionUID=1L;
	
	private final String[] columns=new String[] {"#", "Precursor M/Z", "Charge", "Peptide", "Protein", "Retention Time", "Score"};

	ArrayList<LibraryEntry> entries=new ArrayList<LibraryEntry>();
	
	public void updateEntries(ArrayList<LibraryEntry> newEntries) {
		entries.clear();
		entries.addAll(newEntries);
		fireTableDataChanged();
	}
	
	public LibraryEntry getSelectedRow(int rowIndex) {
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
			case 1: return Double.class;
			case 2: return Byte.class;
			case 3: return String.class;
			case 4: return String.class;
			case 5: return Float.class;
			case 6: return Float.class;
		}
		return Object.class;
	}

	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		LibraryEntry entry=getSelectedRow(rowIndex);
		
		switch (columnIndex) {
			case 0: return rowIndex;
			case 1: return entry.getPrecursorMZ();
			case 2: return entry.getPrecursorCharge();
			case 3: return entry.getPeptideModSeq();
			case 4: return PSMData.accessionsToString(entry.getAccessions());
			case 5: return entry.getRetentionTime()/60f;
			case 6: return entry.getScore();
		}
		return null;
	}
	
	
}
