package edu.washington.gs.maccoss.encyclopedia.gui.dia.interactive;

import java.util.ArrayList;
import java.util.Scanner;
import java.util.StringTokenizer;

import javax.swing.table.AbstractTableModel;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class PeptidePrecursorTableModel extends AbstractTableModel {
	private static final String DELIM = "\t";
	private static final long serialVersionUID=1L;
	private final String[] columns=new String[] {"#", "Peptide", "RT (min)", "Charge", "Passes"};

	private final ArrayList<InteractivePeptidePrecursor> allEntries=new ArrayList<InteractivePeptidePrecursor>();

	public void updateEntries(ArrayList<InteractivePeptidePrecursor> newEntries) {
		allEntries.clear();
		allEntries.addAll(newEntries);
	}
	
	public ArrayList<InteractivePeptidePrecursor> getAllEntries() {
		return allEntries;
	}
	
	public String copy() {
		StringBuilder sb=new StringBuilder(InteractivePeptidePrecursor.getHeader());
		for (int i = 0; i < getRowCount(); i++) {
			sb.append("\n");
			sb.append(getSelectedRow(i).toString());
		}
		return sb.toString();
	}
	
	public void paste(String s) {
		ArrayList<InteractivePeptidePrecursor> entries=new ArrayList<>();
		
		Scanner scanner = new Scanner(s);
		while (scanner.hasNextLine()) {
		  String line = scanner.nextLine();
		  if (line.startsWith("#")) continue; // skip header
		  entries.add(InteractivePeptidePrecursor.parseFromLine(line));
		}
		scanner.close();
		Logger.logLine("Found "+entries.size()+" peptides...");
		
		allEntries.clear();
		allEntries.addAll(entries);
		
		fireTableDataChanged();
	}
	
	public void paste(ArrayList<LibraryEntry> entries) {
		allEntries.clear();
		for (LibraryEntry entry : entries) {
			allEntries.add(new InteractivePeptidePrecursor(entry));
		}
		fireTableDataChanged();
	}
	
	public InteractivePeptidePrecursor getSelectedRow(int rowIndex) {
		return allEntries.get(rowIndex);
	}

	@Override
	public int getRowCount() {
		return allEntries.size();
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
			case 3: return Integer.class;
			case 4: return String.class;
		}
		return Object.class;
	}

	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		InteractivePeptidePrecursor entry=getSelectedRow(rowIndex);
		
		switch (columnIndex) {
			case 0: return rowIndex;
			case 1: return entry.getPeptideModSeq();
			case 2: return entry.getRetentionTimeInSec()/60f;
			case 3: return entry.getPrecursorCharge();
			case 4: 
				if (entry.getIsPassing()==1) {
					return "yes";
				} else if (entry.getIsPassing()==-1) {
					return "no";
				} else {
					return "";
				}
		}
		return null;
	}
	
}
