package edu.washington.gs.maccoss.encyclopedia.gui.massspec;

import java.awt.BorderLayout;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.AbstractTableModel;

import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;

public class FragmentationTable extends JPanel {
	private static final long serialVersionUID=1L;
	private static final DecimalFormat MASS_FORMAT = new DecimalFormat(".#");
	
	public FragmentationTable(LibraryEntry spec, String sequenceKey, SearchParameters params) {
		super(new BorderLayout());
		
		String peptideModSeq=sequenceKey.replaceAll("\\(", "").replaceAll("\\)", "");
		FragmentationModel model=new FragmentationModel(peptideModSeq, params.getAAConstants());
		FragmentIon[] all=model.getPrimaryIonObjects(params.getFragType(), spec.getPrecursorCharge());
		double[] massArray=spec.getMassArray();
		
		ArrayList<FragmentIon> matched=new ArrayList<FragmentIon>();
		for (FragmentIon fragmentIon : all) {
			boolean match=params.getFragmentTolerance().getIndex(massArray, fragmentIon.mass).isPresent();
			if (match) {
				matched.add(fragmentIon);
			}
		}
		
		FragmentIon[] found=matched.toArray(new FragmentIon[matched.size()]);
		
		FragmentationTableModel tableModel=new FragmentationTableModel(all, found, found);

		JTable table=new JTable(tableModel);
		table.setAutoCreateRowSorter(true);

		add(new JScrollPane(table), BorderLayout.CENTER);
	}
	
	public FragmentationTable(PhosphoLocalizationData data, String sequenceKey, SearchParameters params) {
		super(new BorderLayout());
		
		TransitionRefinementData transitionRefinementData=data.getPassingForms().get(sequenceKey);
		
		FragmentIon[] all;
		if (transitionRefinementData!=null) {
			String peptideModSeq=transitionRefinementData.getPeptideModSeq();
			FragmentationModel model=new FragmentationModel(peptideModSeq, params.getAAConstants());
			all=model.getPrimaryIonObjects(params.getFragType(), transitionRefinementData.getPrecursorCharge());
		} else {
			String peptideModSeq=sequenceKey.replaceAll("\\(", "").replaceAll("\\)", "");
			FragmentationModel model=new FragmentationModel(peptideModSeq, params.getAAConstants());
			all=model.getPrimaryIonObjects(params.getFragType(), (byte)3);
		}
		
		FragmentIon[] targets=data.getUniqueTargetFragments().get(sequenceKey);
		Set<FragmentIon> foundSet=data.getUniqueFragmentIons().get(sequenceKey).keySet();
		FragmentIon[] found=foundSet.toArray(new FragmentIon[foundSet.size()]);
		
		FragmentationTableModel tableModel=new FragmentationTableModel(all, targets, found);

		JTable table=new JTable(tableModel);
		table.setAutoCreateRowSorter(true);

		add(new JScrollPane(table), BorderLayout.CENTER);
	}

	class FragmentationTableModel extends AbstractTableModel {
		private static final long serialVersionUID=1L;
		
		IonType[] types;
		FragmentIon[][] ions;
		boolean[][] wasConsidered;
		boolean[][] wasFound;
		
		public FragmentationTableModel(FragmentIon[] all, FragmentIon[] targets, FragmentIon[] found) {
			HashSet<IonType> typeSet=new HashSet<IonType>();
			int maxIonIndex=0;
			for (FragmentIon ion : all) {
				typeSet.add(ion.type);
				if (maxIonIndex<ion.index) {
					maxIonIndex=ion.index;
				}
			}
			types=typeSet.toArray(new IonType[typeSet.size()]);
			Arrays.sort(types);
			
			ions=new FragmentIon[maxIonIndex][];
			wasConsidered=new boolean[maxIonIndex][];
			wasFound=new boolean[maxIonIndex][];
			for (int i=0; i<ions.length; i++) {
				ions[i]=new FragmentIon[types.length];
				wasConsidered[i]=new boolean[types.length];
				wasFound[i]=new boolean[types.length];
			}

			for (FragmentIon ion : all) {
				int typeIndex=Arrays.binarySearch(types, ion.type);
				int ionIndex=ion.index-1;
				ions[ionIndex][typeIndex]=ion;
			}

			for (FragmentIon ion : targets) {
				int typeIndex=Arrays.binarySearch(types, ion.type);
				int ionIndex=ion.index-1;
				wasConsidered[ionIndex][typeIndex]=true;
			}

			for (FragmentIon ion : found) {
				int typeIndex=Arrays.binarySearch(types, ion.type);
				int ionIndex=ion.index-1;
				wasFound[ionIndex][typeIndex]=true;
			}
		}
		
		@Override
		public String getColumnName(int column) {
			if (column==0||column>types.length) return "#";
			return IonType.toString(types[column-1]);
		}
		@Override
		public int getColumnCount() {
			return types.length+2;
		}
		@Override
		public boolean isCellEditable(int rowIndex, int columnIndex) {
			return false;
		}

		@Override
		public int getRowCount() {
			return ions.length;
		}

		@Override
		public Object getValueAt(int rowIndex, int columnIndex) {
			if (columnIndex==0||columnIndex>types.length) return rowIndex+1;
			columnIndex=columnIndex-1;
			
			FragmentIon ion=ions[rowIndex][columnIndex];
			if (ion==null) return "";
			
			String value=MASS_FORMAT.format(ion.mass);
			if (wasFound[rowIndex][columnIndex]) {
				return "<html><b><font color=green>["+value+"]";
			} else if (wasConsidered[rowIndex][columnIndex]) {
				return "<html><b><font color=red>"+value;
			} else {
				return value;
			}
		}
		
	}
}
