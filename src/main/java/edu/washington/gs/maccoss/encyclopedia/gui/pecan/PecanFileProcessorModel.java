package edu.washington.gs.maccoss.encyclopedia.gui.pecan;

import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import javax.swing.table.AbstractTableModel;

import com.google.common.util.concurrent.ThreadFactoryBuilder;

public class PecanFileProcessorModel extends AbstractTableModel {
	private static final long serialVersionUID=1L;
	
	private final String[] columnNames= {"File", "Progress"};
	private final ArrayList<PecanJob> queue=new ArrayList<PecanJob>();
	private final ExecutorService executor;
	
	public PecanFileProcessorModel() {
		ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("Pecan-%d").setDaemon(true).build();
		LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
		executor=new ThreadPoolExecutor(1, 1, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 
	}
	
	public ArrayList<PecanJob> getQueue() {
		return queue;
	}
	
	@Override
	public int getRowCount() {
		return queue.size();
	}
	@Override
	public int getColumnCount() {
		return columnNames.length;
	}
	@Override
	public String getColumnName(int column) {
		return columnNames[column];
	}
	@Override
	public boolean isCellEditable(int rowIndex, int columnIndex) {
		return false;
	}
	
	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		PecanJob job=queue.get(rowIndex);
		switch (columnIndex) {
			case 0:
				return job.getDiaFile().getName();
			case 1:
				return job.getProgressMessage();

			default:
				return null;
		}
	}
	
	public void addJob(PecanJob job) {
		queue.add(job);
		executor.submit(job);
		fireTableDataChanged();
	}
	
	public void fireJobUpdated(PecanJob job) {
		for (int i=0; i<queue.size(); i++) {
			if (job==queue.get(i)) {
				fireTableRowsUpdated(i, i);
			}
		}
	}
}