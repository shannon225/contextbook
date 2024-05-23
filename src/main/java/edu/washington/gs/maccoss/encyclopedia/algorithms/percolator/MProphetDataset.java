package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.util.ArrayList;


public class MProphetDataset {

	private final ArrayList<String> featureNames;
	private final ArrayList<MProphetData> targetPeptideData;
	private final ArrayList<MProphetData> decoyPeptideData;
	public MProphetDataset(ArrayList<String> featureNames, ArrayList<MProphetData> peptideData) {
		this.featureNames = featureNames;
		targetPeptideData=new ArrayList<MProphetData>();
		decoyPeptideData=new ArrayList<MProphetData>();
		for (MProphetData mProphetData : peptideData) {
			if (mProphetData.isDecoy) {
				decoyPeptideData.add(mProphetData);
			} else {
				targetPeptideData.add(mProphetData);
			}
		}
	}
	
	
	public ArrayList<MProphetData> allData() {
		ArrayList<MProphetData> dataset=new ArrayList<>();
		dataset.addAll(targetPeptideData);
		dataset.addAll(decoyPeptideData);
		return dataset;
	}
	
	public ArrayList<float[]> getTargetData() {
		return getDataset(targetPeptideData);
	}
	
	public ArrayList<float[]> getDecoyData() {
		return getDataset(decoyPeptideData);
	}
	
	private ArrayList<float[]> getDataset(ArrayList<MProphetData> dataset) {
		ArrayList<float[]> data=new ArrayList<float[]>();
		for (MProphetData mProphetData : dataset) {
			data.add(mProphetData.data);
		}
		return data;
	}
	
	public ArrayList<MProphetData> getTargetPeptides() {
		return targetPeptideData;
	}
	
	static public class MProphetData implements Comparable<MProphetData> {
		private final String id;
		private final String sequence;
		private final String protein;
		private final float[] data;
		private final boolean isDecoy;
		
		public MProphetData(String id, String sequence, String protein, float[] data, boolean isDecoy) {
			this.id = id;
			this.sequence = sequence;
			this.protein = protein;
			this.data = data;
			this.isDecoy=isDecoy;
		}
		
		@Override
		public int compareTo(MProphetData o) {
			return id.compareTo(o.id);
		}
		@Override
		public int hashCode() {
			return id.hashCode();
		}
		@Override
		public boolean equals(Object obj) {
			if (obj instanceof MProphetData) return compareTo((MProphetData)obj)==0;
			return false;
		}
		
		public float[] getData() {
			return data;
		}
		public String getId() {
			return id;
		}
		public String getProtein() {
			return protein;
		}
		public String getSequence() {
			return sequence;
		}
		public boolean isDecoy() {
			return isDecoy;
		}
	}
}
