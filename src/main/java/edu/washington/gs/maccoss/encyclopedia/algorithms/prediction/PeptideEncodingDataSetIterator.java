package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;

import java.util.ArrayList;
import java.util.List;

import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.DataSetPreProcessor;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;

public class PeptideEncodingDataSetIterator implements DataSetIterator {
    private static final long serialVersionUID = 1L;
    final private int batch;
	private final ArrayList<PeptideEncoding> peptides;

	public PeptideEncodingDataSetIterator(ArrayList<PeptideEncoding> peptides, int batch) {
		this.peptides = peptides;
        this.batch = batch;
	}

    private int curr = 0;
    private DataSetPreProcessor preProcessor;

    private DataSet getDataSet(int i) {
    	return peptides.get(i).encodeDataset();
    }

    @Override
    public synchronized boolean hasNext() {
        return curr < peptides.size();
    }

    @Override
    public synchronized DataSet next() {
        return next(batch);
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    @Override
    public int inputColumns() {
        return getDataSet(0).getFeatures().columns();
    }

    @Override
    public int totalOutcomes() {
        return getDataSet(0).getLabels().columns();
    }

    @Override
    public boolean resetSupported() {
        return true;
    }

    @Override
    public boolean asyncSupported() {
        //Already in memory -> doesn't make sense to prefetch
        return false;
    }

    @Override
    public synchronized void reset() {
        curr = 0;
    }

    @Override
    public int batch() {
        return batch;
    }

    @Override
    public void setPreProcessor(DataSetPreProcessor preProcessor) {
        this.preProcessor = preProcessor;
    }

    @Override
    public List<String> getLabels() {
        return null;
    }

    @Override
    public DataSet next(int num) {
        int end = curr + num;

        List<DataSet> r = new ArrayList<>();
        if (end >= peptides.size())
            end = peptides.size();
        for (; curr < end; curr++) {
            r.add(getDataSet(curr));
        }

        DataSet d = DataSet.merge(r);
        if (preProcessor != null) {
            if (!d.isPreProcessed()) {
                preProcessor.preProcess(d);
                d.markAsPreProcessed();
            }
        }
        return d;
    }

    @Override
    public DataSetPreProcessor getPreProcessor() {
    	return preProcessor;
    }
	
}
