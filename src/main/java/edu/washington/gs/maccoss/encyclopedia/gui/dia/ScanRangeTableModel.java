package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;

import javax.swing.table.AbstractTableModel;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ScanRangeTracker;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.list.array.TFloatArrayList;

public class ScanRangeTableModel extends AbstractTableModel {
	private static final long serialVersionUID=1L;
	
	private final String[] columns=new String[] {"Time", "Start m/z", "Center m/z", "Stop m/z", "Delta m/z"};

	ArrayList<ScoredObject<Range>> ranges=new ArrayList<>();
	
	public void updateEntries(ScanRangeTracker tracker) {
		ranges.clear();
		HashMap<Range, TFloatArrayList> map=tracker.getStripeRTsInSecs();
		for (Entry<Range, TFloatArrayList> entry : map.entrySet()) {
			Range r=entry.getKey();
			for (float f : entry.getValue().toArray()) {
				ranges.add(new ScoredObject<Range>(f, r));
			}
		}
		Collections.sort(ranges);
		fireTableDataChanged();
	}
	
	public ScoredObject<Range> getSelectedRow(int rowIndex) {
		return ranges.get(rowIndex);
	}
	
	public float getCenter(int rowIndex) {
		ScoredObject<Range> entry=getSelectedRow(rowIndex);
		return entry.y.getMiddle();
	}

	@Override
	public int getRowCount() {
		return ranges.size();
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
			case 0: return Float.class;
			case 1: return Float.class;
			case 2: return Float.class;
			case 3: return Float.class;
			case 4: return Float.class;
		}
		return Object.class;
	}

	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		ScoredObject<Range> entry=getSelectedRow(rowIndex);
		
		switch (columnIndex) {
			case 0: return entry.x;
			case 1: return entry.y.getStart();
			case 2: return entry.y.getMiddle();
			case 3: return entry.y.getStop();
			case 4: return entry.y.getRange();
		}
		return null;
	}
	
	
}
