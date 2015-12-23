package edu.washington.gs.maccoss.encyclopedia.utils.threading;

import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.utils.io.OutputMessage;

public class ExternalExecutorTest {
	public static void main(String[] args) throws Exception {
		ExternalExecutor e=new ExternalExecutor(new String[] {"cat", "/Users/searleb/Documents/data/dbs/UP000005640_9606.fasta"});
		BlockingQueue<OutputMessage> result=e.start();

		int count=0;
		while (!e.isFinished()||!result.isEmpty()) {
			OutputMessage data=result.take();
			if (data.message.startsWith(">")) {
				count++;
				if (count%1000==0) System.out.println(count+", "+result.size()+", "+e.isFinished());
			}
		}
		System.out.println(count+", "+result.size());
		System.out.println("FINISHED!");
	}
}
