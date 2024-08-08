package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;

import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.api.buffer.DataType;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class PeptideEncoding {
    public static final DataType DEFAULT_DATA_TYPE = DataType.FLOAT;
	public static final int MAX_CHARGE=6;
    public static final int MAX_PEPTIDE_LENGTH=32; // number of termini + 30
    
    public static final int ENCODED_INPUT_SIZE=MAX_PEPTIDE_LENGTH*EncodedAminoAcid.MAX_ENCODING_LENGTH+MAX_CHARGE;
    public static final int ENCODED_OUTPUT_SIZE=MAX_PEPTIDE_LENGTH*4+2;
    
	private final String peptideModSeq;
	private final byte charge;
	private final float rtInSec;
	private final float ims;
	private final float[] bp1=new float[MAX_PEPTIDE_LENGTH];
	private final float[] bp2=new float[MAX_PEPTIDE_LENGTH];
	private final float[] yp1=new float[MAX_PEPTIDE_LENGTH];
	private final float[] yp2=new float[MAX_PEPTIDE_LENGTH];
	
	public DataSet encodeDataset(SearchParameters parameters) {
		EncodedAminoAcid[] aas=EncodedAminoAcid.getAAs(peptideModSeq, parameters.getAAConstants());
		//System.out.println("Input: "+encodeInput(aas).shapeInfoToString());
		//System.out.println("Output: "+encodeResult().shapeInfoToString());
		//System.out.println("Input Mask: "+encodeInputMask(aas).shapeInfoToString());
		//System.out.println("Output Mask: "+encodeResultMask(aas).shapeInfoToString());
		
		//Input: Rank: 2, DataType: FLOAT, Offset: 0, Order: c, Shape: [1,1446],  Stride: [1,1]
		//Output: Rank: 2, DataType: FLOAT, Offset: 0, Order: c, Shape: [1,130],  Stride: [1,1]
		//Input Mask: Rank: 2, DataType: FLOAT, Offset: 0, Order: c, Shape: [1,1446],  Stride: [1,1]
		//Output Mask: Rank: 2, DataType: FLOAT, Offset: 0, Order: c, Shape: [1,130],  Stride: [1,1]
						
		// TODO why does the input mask need to be absent? 
    	return new DataSet(encodeInput(aas), encodeResult(), null, encodeResultMask(aas));
    	//return new DataSet(encodeInput(aas), encodeResult());
	}
	
	public INDArray encodeInputMask(EncodedAminoAcid[] aas) {
    	float[] mask=new float[ENCODED_INPUT_SIZE];
        int start=aas[0].isNTerm()?0:1;
        for (int i = start; i < aas.length; i++) {
        	for (int j = 0; j < MAX_PEPTIDE_LENGTH; j++) {
            	int index=i*MAX_PEPTIDE_LENGTH+j;
            	mask[index]=1.0f;
        	}
        }
        for (int i = MAX_PEPTIDE_LENGTH*EncodedAminoAcid.MAX_ENCODING_LENGTH; i < ENCODED_INPUT_SIZE; i++) {
        	mask[i]=1.0f;
		}
        
        return Nd4j.create(mask, new long[] {1, ENCODED_INPUT_SIZE}, DEFAULT_DATA_TYPE);
	}

	public INDArray encodeInput(EncodedAminoAcid[] aas) {
    	float[] data=new float[ENCODED_INPUT_SIZE];

        int start=aas[0].isNTerm()?0:1;
        for (int i = start; i < aas.length; i++) {
        	int index=i*MAX_PEPTIDE_LENGTH+aas[i].getIndex();
        	data[index]=1.0f;
        }

        // set "zero" for empty amino acids
        for (int i = 0; i < start; i++) {
        	int index=i*MAX_PEPTIDE_LENGTH;
        	data[index]=1.0f;
        }

        // set "zero" for empty amino acids
        for (int i = aas.length; i < MAX_PEPTIDE_LENGTH; i++) {
        	int index=i*MAX_PEPTIDE_LENGTH;
        	data[index]=1.0f;
        }
        
        // set charge
        data[MAX_PEPTIDE_LENGTH*EncodedAminoAcid.MAX_ENCODING_LENGTH+charge-1]=1.0f;
    	
        return Nd4j.create(data, new long[] {1, ENCODED_INPUT_SIZE}, DEFAULT_DATA_TYPE);
	}
	
	public INDArray encodeResultMask(EncodedAminoAcid[] aas) {
		float[] mask=new float[MAX_PEPTIDE_LENGTH];
        for (int i = 0; i < aas.length; i++) {
        	mask[i]=1.0f;
        }

		float[] data=General.concatenate(mask, mask, mask, mask, new float[] {1.0f, 1.0f});
        return Nd4j.create(data, new long[] {1, ENCODED_OUTPUT_SIZE}, DEFAULT_DATA_TYPE);
	}
	
	public INDArray encodeResult() {
		float[] data=General.concatenate(bp1, bp2, yp1, yp2, new float[] {rtInSec, ims});
		return Nd4j.create(data, new long[] {1, ENCODED_OUTPUT_SIZE}, DEFAULT_DATA_TYPE);
	}
	
	public PeptideEncoding(LibraryEntry entry, float rtInSec, SearchParameters parameters) {
		this.peptideModSeq=entry.getPeptideModSeq();
		this.charge=entry.getPrecursorCharge();
		this.rtInSec=rtInSec;
		this.ims=entry.getIonMobility().isEmpty()?-1.0f:entry.getIonMobility().get();
		
		double[] massArray = entry.getMassArray();
		float[] intensityArray=entry.getIntensityArray();
		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		for (FragmentIon fragmentIon : model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), true)) {
			byte fragcharge=IonType.getCharge(fragmentIon.getType());
			if (fragcharge>3||fragcharge<1) {
				continue;
			}
			IonType type=IonType.getCanonicalIonType(fragmentIon.getType());
			float[] array=null;
			if (type==IonType.b) {
				switch (fragcharge) {
					case 1: array=bp1; break;
					case 2: array=bp2; break;
					default: break;
				}
			} else if (type==IonType.y) {
				switch (fragcharge) {
					case 1: array=yp1; break;
					case 2: array=yp2; break;
					default: break;
				}
			}
			if (array==null) {
				continue;
			}
			
			float intensity=parameters.getFragmentTolerance().getMaxIntensity(massArray, intensityArray, fragmentIon.getMass());
			int fragindex=fragmentIon.getIndex();
			
			array[fragindex]=intensity;
		}
	}
}
