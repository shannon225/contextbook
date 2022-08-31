package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import javax.swing.table.AbstractTableModel;

import com.google.common.util.concurrent.ThreadFactoryBuilder;

import edu.washington.gs.maccoss.encyclopedia.jobs.JobProcessor;
import edu.washington.gs.maccoss.encyclopedia.jobs.WorkerJob;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class JobProcessorTableModel extends AbstractTableModel implements JobProcessor {
	private static final long serialVersionUID=1L;
	
	private final String[] columnNames= {"File", "Progress"};
	private final ArrayList<SwingWorkerJob> queue=new ArrayList<SwingWorkerJob>();
	private final ExecutorService executor;
	
	public JobProcessorTableModel() {
		ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("EncyclopeDIA-%d").setDaemon(true).build();
		LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
		executor=new ThreadPoolExecutor(1, 1, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 
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
		SwingJob job=queue.get(rowIndex);
		switch (columnIndex) {
			case 0:
				return job.getJobTitle();
			case 1:
				return job.getProgressMessage();

			default:
				return null;
		}
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.gui.pecan.JobProcessor#getQueue()
	 */
	@Override
	public ArrayList<WorkerJob> getQueue() {
		ArrayList<WorkerJob> jobs=new ArrayList<WorkerJob>();
		for (SwingWorkerJob job : queue) {
			jobs.add(job.getJob());
		}
		return jobs;
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.gui.pecan.JobProcessor#addJob(edu.washington.gs.maccoss.encyclopedia.gui.pecan.PecanJob)
	 */
	@Override
	public void addJob(WorkerJob job) {
		Logger.logLine("Adding new job to queue: "+job.getJobTitle());
		SwingWorkerJob swingJob=new SwingWorkerJob(this, job);
		queue.add(swingJob);
		executor.submit(swingJob);
		fireTableDataChanged();
	}
	
	public void clearJobs() {
		queue.clear();
		fireTableDataChanged();
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.gui.pecan.JobProcessor#fireJobUpdated(edu.washington.gs.maccoss.encyclopedia.gui.pecan.PecanJob)
	 */
	@Override
	public void fireJobUpdated(WorkerJob job) {
		for (int i=0; i<queue.size(); i++) {
			if (job==queue.get(i).getJob()) {
				fireTableRowsUpdated(i, i);
			}
		}
	}
}