package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Optional;

import org.nd4j.linalg.api.buffer.DataType;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.MultiDataSet;
import org.nd4j.linalg.factory.Nd4j;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class PeptideEncoding implements MultiDataSetEncoding {
    private static final float MINIMUM_INTENSITY_THRESHOLD_FOR_TRAINING = 0.05f;
	public static final DataType DEFAULT_DATA_TYPE = DataType.FLOAT;
	public static final int MAX_CHARGE=6;
    public static final int MAX_PEPTIDE_LENGTH=32; // number of termini + 30
    
    public static final int ENCODED_INPUT_PEPTIDE_SIZE=MAX_PEPTIDE_LENGTH*AminoAcidEncoding.MAX_ENCODING_LENGTH;
    public static final int ENCODED_INPUT_CHARGE_SIZE=MAX_CHARGE;
    
    public static final int ENCODED_OUTPUT_FRAGMENT_SIZE=MAX_PEPTIDE_LENGTH*4;
    public static final int ENCODED_OUTPUT_IMS_SIZE=1;
    public static final int ENCODED_OUTPUT_RT_SIZE=1;
    
	private final AminoAcidEncoding[] aas;
	private final byte charge;
	private final float rtInSec;
	private final float ims;
	private final float[] bp1;
	private final float[] bp2;
	private final float[] yp1;
	private final float[] yp2;
	
	public MultiDataSet encodeDataset() {
		//System.out.println("Input: "+encodeInput(aas).shapeInfoToString());
		//System.out.println("Output: "+encodeResult().shapeInfoToString());
		//System.out.println("Input Mask: "+encodeInputMask(aas).shapeInfoToString());
		//System.out.println("Output Mask: "+encodeResultMask(aas).shapeInfoToString());
		
		//Input: Rank: 2, DataType: FLOAT, Offset: 0, Order: c, Shape: [1,1446],  Stride: [1,1]
		//Output: Rank: 2, DataType: FLOAT, Offset: 0, Order: c, Shape: [1,130],  Stride: [1,1]
		//Input Mask: Rank: 2, DataType: FLOAT, Offset: 0, Order: c, Shape: [1,1446],  Stride: [1,1]
		//Output Mask: Rank: 2, DataType: FLOAT, Offset: 0, Order: c, Shape: [1,130],  Stride: [1,1]
						
		// TODO why does the input mask need to be absent? 
		return new MultiDataSet(encodeInput(), encodeResult(), null, encodeResultMask());
	}

	public INDArray[] encodeInput() {
		return encodeInput(aas, charge);
	}
	
	public static INDArray[] encodeInput(String peptideModSeq, byte charge, AminoAcidConstants constants) {
		AminoAcidEncoding[] aas=AminoAcidEncoding.getAAs(peptideModSeq, constants);
		return encodeInput(aas, charge);
	}

	public static INDArray[] encodeInput(AminoAcidEncoding[] aas, byte charge) {
    	float[] data=new float[ENCODED_INPUT_PEPTIDE_SIZE];

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
        float[] chargeArray=new float[ENCODED_INPUT_CHARGE_SIZE];
        Arrays.fill(chargeArray, -1.0f);
        chargeArray[charge-1]=1.0f;
        
        return new INDArray[] {
        		Nd4j.create(data, new long[] {1, AminoAcidEncoding.MAX_ENCODING_LENGTH*MAX_PEPTIDE_LENGTH}, DEFAULT_DATA_TYPE), 
        		Nd4j.create(chargeArray, new long[] {1, ENCODED_INPUT_CHARGE_SIZE}, DEFAULT_DATA_TYPE)
        };
	}
	
	public INDArray[] encodeResultMask() {
		float[] mask=new float[MAX_PEPTIDE_LENGTH];
        for (int i = 0; i < aas.length; i++) {
        	mask[i]=1.0f;
        }

		float[] data=General.concatenate(mask, mask, mask, mask);

		return new INDArray[] {
				Nd4j.create(data, new long[] {1, ENCODED_OUTPUT_FRAGMENT_SIZE}, DEFAULT_DATA_TYPE),
				Nd4j.create(new float[] {1.0f}, new long[] {1, ENCODED_OUTPUT_RT_SIZE}, DEFAULT_DATA_TYPE),
				Nd4j.create(new float[] {1.0f}, new long[] {1, ENCODED_OUTPUT_IMS_SIZE}, DEFAULT_DATA_TYPE)
		};
	}
	
	public INDArray[] encodeResult() {
		float[] data=new float[ENCODED_OUTPUT_FRAGMENT_SIZE];
		
		System.arraycopy(bp1, 0, data, 0, bp1.length);
		System.arraycopy(bp2, 0, data, MAX_PEPTIDE_LENGTH, bp2.length);
		System.arraycopy(yp1, 0, data, MAX_PEPTIDE_LENGTH*2, yp1.length);
		System.arraycopy(yp2, 0, data, MAX_PEPTIDE_LENGTH*3, yp2.length);
		
		return new INDArray[] {
				Nd4j.create(data, new long[] {1, ENCODED_OUTPUT_FRAGMENT_SIZE}, DEFAULT_DATA_TYPE),
				Nd4j.create(new float[] {rtInSec}, new long[] {1, ENCODED_OUTPUT_RT_SIZE}, DEFAULT_DATA_TYPE),
				Nd4j.create(new float[] {ims}, new long[] {1, ENCODED_OUTPUT_IMS_SIZE}, DEFAULT_DATA_TYPE)
		};
	}
	
	public PeptideEncoding(LibraryEntry entry, float rtInSec, SearchParameters parameters) {
		this.aas=AminoAcidEncoding.getAAs(entry.getPeptideModSeq(), parameters.getAAConstants());
		this.charge=entry.getPrecursorCharge();
		this.rtInSec=rtInSec;
		this.ims=entry.getIonMobility().isEmpty()?-1.0f:entry.getIonMobility().get();
		
		double[] massArray = entry.getMassArray();
		float[] intensityArray=entry.getIntensityArray();
		float intensityThreshold=General.max(intensityArray)*MINIMUM_INTENSITY_THRESHOLD_FOR_TRAINING;
		
		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		bp1=new float[aas.length];
		bp2=new float[aas.length];
		yp1=new float[aas.length];
		yp2=new float[aas.length];
		
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
			int fragindex=fragmentIon.getIndex()-1;
			
			if (intensity>intensityThreshold) {
				array[fragindex]=intensity;
			}
		}
	}
	
	public float score(INDArray[] output) {
		INDArray[] expected=encodeResult();
		
		throw new EncyclopediaException("Not IMPLEMENTED YET!!!!"); //FIXME
		//return Correlation.getSpectralContrastAngle(expected.toFloatVector(), output.toFloatVector());
	}
	
	public LibraryEntry encodingToEntry(HashSet<String> accessions, AminoAcidConstants constants) {
		return outputToEntry(AminoAcidEncoding.getPeptideModSeq(aas, constants), charge, accessions, encodeResult(), constants);
	}
	
	public LibraryEntry outputToEntry(HashSet<String> accessions, INDArray[] output, AminoAcidConstants constants) {
		return outputToEntry(AminoAcidEncoding.getPeptideModSeq(aas, constants), charge, accessions, output, constants);
	}
	
	public static LibraryEntry outputToEntry(String peptideModSeq, byte charge, HashSet<String> accessions, INDArray[] output, AminoAcidConstants constants) {
		float[] data=output[0].toFloatVector();
		float rtInSec=output[1].getFloat(0);
		float ionMobility=output[2].getFloat(0);
		
		return outputToEntry(peptideModSeq, data, rtInSec, ionMobility, charge, accessions, constants);
	}

	public static LibraryEntry outputToEntry(String peptideModSeq, float[] data, float rtInSec, float ionMobility,
			byte charge, HashSet<String> accessions, AminoAcidConstants constants) {
		FragmentationModel model=PeptideUtils.getPeptideModel(peptideModSeq, constants);
		FragmentIon[] bIons=model.getBIons(false);
		FragmentIon[] bp2Ions=FragmentationModel.getPlus2s(bIons);
		FragmentIon[] yIons=model.getYIons(false);
		FragmentIon[] yp2Ions=FragmentationModel.getPlus2s(yIons);
		
		ArrayList<Peak> peaks=new ArrayList<Peak>();
		
		int index=0;
		for (int i = 0; i < bIons.length; i++) {
			if (data[index]>0) {
				Peak p=new Peak(bIons[i].getMass(), data[index]);
				peaks.add(p);
			}
			index++;
		}
		index=MAX_PEPTIDE_LENGTH;
		for (int i = 0; i < bp2Ions.length; i++) {
			if (data[index]>0) {
				Peak p=new Peak(bp2Ions[i].getMass(), data[index]);
				peaks.add(p);
			}
			index++;
		}
		index=MAX_PEPTIDE_LENGTH*2;
		for (int i = 0; i < yIons.length; i++) {
			if (data[index]>0) {
				Peak p=new Peak(yIons[i].getMass(), data[index]);
				peaks.add(p);
			}
			index++;
		}
		index=MAX_PEPTIDE_LENGTH*3;
		for (int i = 0; i < yp2Ions.length; i++) {
			if (data[index]>0) {
				Peak p=new Peak(yp2Ions[i].getMass(), data[index]);
				peaks.add(p);
			}
			index++;
		}
		
		Collections.sort(peaks);
		Triplet<double[], float[], Optional<float[]>> triplet=Peak.toArrays(peaks);
		double[] massArray=triplet.x;
		float[] intensityArray=triplet.y;
		
		
		return new LibraryEntry(peptideModSeq, accessions, model.getChargedMass(charge), charge, peptideModSeq, 0, rtInSec, 0.0f, massArray, intensityArray, Optional.of(ionMobility), constants);
	}
}
