package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;
import java.util.concurrent.BlockingQueue;

import junit.framework.TestCase;

public class PercolatorExecutorTest extends TestCase {

	public void testPercolatorExecutor() throws Exception {
		PercolatorExecutor e=new PercolatorExecutor(new File("/Users/searleb/Documents/projects/pecan/v0.9.7/onePep.530.49.td.feature"));
		BlockingQueue<String> result=e.start();

		int count=0;
		while (!e.isFinished()||!result.isEmpty()) {
			String data=result.take();
			if (data.startsWith(">")) {
				count++;
				if (count%1000==0) System.out.println(count+", "+result.size()+", "+e.isFinished());
			}
		}
		System.out.println(count+", "+result.size());
		System.out.println("FINISHED!");
	}
}
