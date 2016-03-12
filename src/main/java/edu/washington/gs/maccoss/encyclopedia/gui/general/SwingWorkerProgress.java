package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.awt.BorderLayout;
import java.awt.Frame;
import java.util.concurrent.ExecutionException;

import javax.swing.BorderFactory;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.SwingWorker;

import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;

public abstract class SwingWorkerProgress<T> {
	public static void main(String[] args) {
		SwingWorkerProgress<Nothing> worker=new SwingWorkerProgress<Nothing>(null, "Working Title...", "Working text...") {
			@Override
			protected Nothing doInBackgroundForReal() throws Exception {
				System.out.println("SHOULD HAVE STARTED");
				Thread.sleep(1000);
				System.out.println("1...");
				Thread.sleep(1000);
				System.out.println("2...");
				Thread.sleep(1000);
				System.out.println("3...");
				Thread.sleep(1000);
				return Nothing.NOTHING;
			}
			@Override
			protected void doneForReal(Nothing t) {
				System.out.println("SHOULD BE DONE");
			}
		};
		
		worker.execute();
	}
	final JDialog dlgProgress;
	final SwingWorker<T, Void> sw;
	public SwingWorkerProgress(Frame frame, String title, String text) {
		dlgProgress=new JDialog(frame, title, true);

		JProgressBar pbProgress=new JProgressBar(0, 100);
		pbProgress.setIndeterminate(true);

		JPanel panel=new JPanel(new BorderLayout());
		panel.add(BorderLayout.CENTER, pbProgress);
		panel.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10), BorderFactory.createTitledBorder(text)));
		
		dlgProgress.add(BorderLayout.CENTER, panel);
		dlgProgress.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
		dlgProgress.setSize(300, 90);

		sw=new SwingWorker<T, Void>() {
			@Override
			protected T doInBackground() throws Exception {
				return doInBackgroundForReal();
			}

			@Override
			protected void done() {
				try {
					doneForReal(get());
				} catch (InterruptedException e) {
					System.err.println("Encountered Fatal Error reading library!");
					e.printStackTrace();
				} catch (ExecutionException e) {
					System.err.println("Encountered Fatal Error reading library!");
					e.printStackTrace();
				}
				dlgProgress.setVisible(false);
				dlgProgress.dispose();
			}
		};

	}
	
	public void execute() {
		sw.execute();
		dlgProgress.setVisible(true);
	}

	protected abstract T doInBackgroundForReal() throws Exception;

	protected abstract void doneForReal(T t);
}
