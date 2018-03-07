package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.util.ArrayList;
import java.util.Map.Entry;
import java.util.TreeMap;

import javax.swing.table.AbstractTableModel;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LocalizedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;

public class LocalizedLibraryEntryTableModel extends AbstractTableModel {
	private static final long serialVersionUID=1L;
	
	public static final int deltaRTColumnIndex=5;
	private final String[] columns=new String[] {"#", "Peptide", "Number of Mods", "Number of Forms", "Number of Sites", "Delta RT (sec)", "Protein"};

	ArrayList<Pair<String, ArrayList<LocalizedLibraryEntry>>> entries=new ArrayList<Pair<String, ArrayList<LocalizedLibraryEntry>>>();
	
	public void updateEntries(ArrayList<LocalizedLibraryEntry> newEntries) {
		entries.clear();
		TreeMap<String, ArrayList<LocalizedLibraryEntry>> entryMap=new TreeMap<>();
		for (LocalizedLibraryEntry entry : newEntries) {
			if (entry.getNumberOfModifiableResidues()==1||entry.getLocalizationIons().length>0) {
				String key = entry.getPeptideSeq()+"+"+entry.getNumberOfModifications()+"m";
				ArrayList<LocalizedLibraryEntry> list=entryMap.get(key);
				if (list==null) {
					list=new ArrayList<>();
					entryMap.put(key, list);
				}
				list.add(entry);
			}
		}
		
		for (Entry<String, ArrayList<LocalizedLibraryEntry>> entry : entryMap.entrySet()) {
			entries.add(new Pair<String, ArrayList<LocalizedLibraryEntry>>(entry.getKey(), entry.getValue()));
		}
		
		fireTableDataChanged();
	}
	
	public ArrayList<LocalizedLibraryEntry> getSelectedRow(int rowIndex) {
		return entries.get(rowIndex).y;
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
			case 2: return Integer.class;
			case 3: return Integer.class;
			case 4: return Integer.class;
			case 5: return Float.class;
			case 6: return String.class;
		}
		return Object.class;
	}
	
	public boolean flagRow(int rowIndex) {
		ArrayList<LocalizedLibraryEntry> value=entries.get(rowIndex).y;
		if (value.size()<=1) return false;
		
		float deltaRT = getDeltaRT(value);
		return deltaRT<5;
	}

	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		String key=entries.get(rowIndex).x;
		ArrayList<LocalizedLibraryEntry> value=entries.get(rowIndex).y;
		float deltaRT = getDeltaRT(value);
		
		switch (columnIndex) {
			case 0: return rowIndex;
			case 1: return key;
			case 2: return value.get(0).getNumberOfModifications();
			case 3: return value.size();
			case 4: return value.get(0).getNumberOfModifiableResidues();
			case 5: return deltaRT;
			case 6: return PSMData.accessionsToString(value.get(0).getAccessions());
		}
		return null;
	}

	private float getDeltaRT(ArrayList<LocalizedLibraryEntry> value) {
		float minRT=Float.MAX_VALUE;
		float maxRT=-Float.MAX_VALUE;
		for (LocalizedLibraryEntry entry : value) {
			if (entry.getRetentionTime()>maxRT) maxRT=entry.getRetentionTime();
			if (entry.getRetentionTime()<minRT) minRT=entry.getRetentionTime();
		}
		float deltaRT=maxRT-minRT;
		return deltaRT;
	}
	
	
}
