package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.util.ArrayList;

import javax.swing.table.AbstractTableModel;

import edu.washington.gs.maccoss.encyclopedia.filewriters.LibraryReportExtractor.PeptideReportData;

public class MultiPeptideResultsTableModel extends AbstractTableModel {
	private static final long serialVersionUID=1L;
	
	private final String[] columns=new String[] {"#", "Peptide", "Protein", "Quant Ions", "Average RT"};

	ArrayList<PeptideReportData> allEntries=new ArrayList<PeptideReportData>();
	ArrayList<PeptideReportData> entries=new ArrayList<PeptideReportData>();
	private int minimumNumberOfTransitions=1;
	
	public void updateEntries(ArrayList<PeptideReportData> newEntries) {
		allEntries.clear();
		allEntries.addAll(newEntries);
		
		filterTable(minimumNumberOfTransitions);
	}
	
	public void filterTable(int minimumNumberOfTransitions) {
		entries.clear();
		this.minimumNumberOfTransitions=minimumNumberOfTransitions;
		for (PeptideReportData entry : allEntries) {
			if (minimumNumberOfTransitions<=entry.getTargetFragmentMzs().length) {
				entries.add(entry);
			}
		}
		
		fireTableDataChanged();
	}
	
	public PeptideReportData getSelectedRow(int rowIndex) {
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
			case 2: return String.class;
			case 3: return Integer.class;
			case 4: return Float.class;
		}
		return Object.class;
	}

	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		PeptideReportData entry=getSelectedRow(rowIndex);
		
		switch (columnIndex) {
			case 0: return rowIndex;
			case 1: return entry.getPeptideModSeq();
			case 2: return entry.getAccessions();
			case 3: return entry.getTargetFragmentMzs().length;
			case 4: return entry.getAverageRetentionTime()/60f;
		}
		return null;
	}
	
	
}
