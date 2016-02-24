package edu.washington.gs.maccoss.encyclopedia.utils.math.randomforest;

import gnu.trove.map.hash.TIntIntHashMap;

import java.util.ArrayList;

public class Forest {
	ArrayList<TreeNode> trees=new ArrayList<TreeNode>();

	public Forest() {
	}

	public void add(TreeNode t) {
		trees.add(t);
	}

	public boolean predictValue(Record record) {
		float probability=predictProbability(record);
		return probability>0.5f;
	}

	public float predictProbability(Record record) {
		int sum=0;
		for (TreeNode tree : trees) {
			if (tree.predictValue(record)) {
				sum++;
			}
		}
		float probability=sum/(float) trees.size();
		return probability;
	}

	public void addDecisionIndicies(TIntIntHashMap map) {
		for (TreeNode node : trees) {
			node.addDecisionIndicies(map);
		}
	}
}
