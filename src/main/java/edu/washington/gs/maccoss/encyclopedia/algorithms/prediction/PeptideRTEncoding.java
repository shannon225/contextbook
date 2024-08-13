package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.factory.Nd4j;

public class PeptideRTEncoding implements DataSetEncoding {
    public static final int MAXIMUM_PEPTIDE_LENGTH = 50;
	private final AminoAcidEncoding[] aas;
	private final float hirt;

	public PeptideRTEncoding(AminoAcidEncoding[] aas, float hirt) {
		this.aas = aas;
		this.hirt=hirt;
	}
	
	@Override
	public DataSet encodeDataset() {
        // One-hot encode the sequences
		INDArray input= AminoAcidEncoding.encode(MAXIMUM_PEPTIDE_LENGTH, aas);
        INDArray output = Nd4j.create(new double[][]{new double[] {hirt}});
        return new DataSet(input, output);
	}
	
	public float getHiRT() {
		return hirt;
	}
}
