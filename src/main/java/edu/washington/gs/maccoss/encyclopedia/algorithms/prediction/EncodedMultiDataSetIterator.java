package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;

import java.util.ArrayList;
import java.util.List;

import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.MultiDataSetPreProcessor;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

public class EncodedMultiDataSetIterator implements MultiDataSetIterator {
    private static final long serialVersionUID = 1L;
    final private int batch;
	private final ArrayList<? extends MultiDataSetEncoding> list;

	private int curr = 0;
    private MultiDataSetPreProcessor preProcessor;

	public EncodedMultiDataSetIterator(ArrayList<? extends MultiDataSetEncoding> list, int batch) {
		this.list = list;
        this.batch = batch;
	}
	

    @Override
    public MultiDataSet next(int num) {
        int end = curr + num;

        List<MultiDataSet> r = new ArrayList<>();
        if (end >= list.size()) {
            end = list.size();
        }
        for (; curr < end; curr++) {
            r.add(list.get(curr).encodeDataset());
        }

        MultiDataSet d = org.nd4j.linalg.dataset.MultiDataSet.merge(r);
        if (preProcessor != null) {
            preProcessor.preProcess(d);
        }
        return d;
    }

    @Override
    public void setPreProcessor(MultiDataSetPreProcessor preProcessor) {
        this.preProcessor = preProcessor;
    }

    @Override
    public MultiDataSetPreProcessor getPreProcessor() {
        return this.preProcessor;
    }

    @Override
    public boolean resetSupported() {
        return true;
    }

    @Override
    public boolean asyncSupported() {
        return false;
    }

    @Override
    public void reset() {
        curr = 0;
    }

    @Override
    public boolean hasNext() {
        return curr < list.size();
    }

    @Override
    public MultiDataSet next() {
        return next(batch);
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }
	
}
