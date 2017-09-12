package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;

import javax.swing.table.AbstractTableModel;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Map;
import java.util.Map.Entry;

public class SampleTableModel extends AbstractTableModel {
	private static final long serialVersionUID=1L;
	
	private final String[] columns=new String[] {"#", "Sample", "Min", "Max", "Range", "Cycle", "TIC"};

	ArrayList<Pair<StripeFileInterface, float[]>> entries=new ArrayList<Pair<StripeFileInterface,float[]>>();
	
	public void updateEntries(ArrayList<StripeFileInterface> newEntries) {
		entries.clear();
		for (StripeFileInterface file : newEntries) {
			Map<Range, Float> ranges=file.getRanges();
			
			TFloatArrayList cycleTimes=new TFloatArrayList();
			TFloatArrayList rangeWindows=new TFloatArrayList();
			float min=Float.MAX_VALUE;
			float max=0.0f;
			for (Entry<Range, Float> entry : ranges.entrySet()) {
				cycleTimes.add(entry.getValue());
				Range key=entry.getKey();
				if (key.getStart()<min) {
					min=key.getStart();
				}
				if (key.getStop()>max) {
					max=key.getStop();
				}
				rangeWindows.add(key.getRange());
			}
			float[] sampleData=new float[] {min, max, General.mean(rangeWindows.toArray()), General.mean(cycleTimes.toArray())};
			
			entries.add(new Pair<StripeFileInterface, float[]>(file, sampleData));
		}
		fireTableDataChanged();
	}
	
	public StripeFileInterface getSelectedRow(int rowIndex) {
		return entries.get(rowIndex).x;
	}
	
	public ArrayList<StripeFileInterface> getRows() {
		ArrayList<StripeFileInterface> files=new ArrayList<StripeFileInterface>();
		for (Pair<StripeFileInterface, float[]> pair : entries) {
			files.add(pair.x);
		}
		return files;
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
			case 3: return Float.class;
			case 4: return Float.class;
			case 5: return Float.class;
			case 6: return Float.class;
		}
		return Object.class;
	}

	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		Pair<StripeFileInterface, float[]> pair=entries.get(rowIndex);
		
		switch (columnIndex) {
			case 0: return rowIndex;
			case 1: return pair.x.getOriginalFileName();
			case 2: return pair.y[0];
			case 3: return pair.y[1];
			case 4: return pair.y[2];
			case 5: return pair.y[3];
			case 6: {
				try {
					return pair.x.getTIC();
				} catch (IOException ioe) {
					Logger.errorLine("Error reading raw file.");
					Logger.errorException(ioe);
				} catch (SQLException sqle) {
					Logger.errorLine("Error reading raw file.");
					Logger.errorException(sqle);
				}
			}
		}
		return null;
	}
	
	
}
