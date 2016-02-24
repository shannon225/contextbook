package edu.washington.gs.maccoss.encyclopedia.utils.math.randomforest;

import gnu.trove.map.hash.TIntIntHashMap;

public class TreeNode implements Comparable<TreeNode> {
	private final TreeNode parent;
	private TreeNode left=null;
	private TreeNode right=null;

	private final DecisionTreeSplitPoint decision;
	private final float probability;

	public TreeNode(TreeNode parent, DecisionTreeSplitPoint decision, float probability) {
		this.parent=parent;
		this.decision=decision;
		this.probability=probability;
	}
	
	@Override
	public int compareTo(TreeNode o) {
		if (o==null) return 1;
		int c=decision.compareTo(o.decision);
		if (c!=0) return 0;
		c=Float.compare(probability, o.probability);
		if (c!=0) return 0;
		c=compareTo(left, o.left);
		if (c!=0) return 0;
		c=compareTo(right, o.right);
		return c;
	}
	
	public static int compareTo(TreeNode o1, TreeNode o2) {
		if (o1==null&&o2==null) return 0;
		if (o1==null) return -1;
		if (o2==null) return 1;
		return o1.compareTo(o2);
	}

	public void addDecisionIndicies(TIntIntHashMap map) {
		if (decision!=null) {
			map.adjustOrPutValue(decision.index, 1, 1);
			if (left!=null)
				left.addDecisionIndicies(map);
			if (right!=null)
				right.addDecisionIndicies(map);
		}
	}

	public TreeNode getParent() {
		return parent;
	}

	public TreeNode getLeft() {
		return left;
	}

	public TreeNode getRight() {
		return right;
	}

	public TreeNode getRandomLeaf() {
		if (Math.random()>0.5f) {
			if (left!=null) {
				return left.getRandomLeaf();
			}
		} else {
			if (right!=null) {
				return right.getRandomLeaf();
			}
		}
		return this;
	}

	public boolean predictValue(Record record) {
		if (isLeaf())
			return probability>0.5f;
		if (decision.isRight(record)) {
			if (right==null)
				return probability>0.5f;
			return right.predictValue(record);
		} else {
			if (left==null)
				return probability>0.5f;
			return left.predictValue(record);
		}
	}

	public boolean isLeaf() {
		return decision==null;
	}

	public void setLeft(TreeNode left) {
		this.left=left;
	}

	public void setRight(TreeNode right) {
		this.right=right;
	}

	@Override
	public String toString() {
		return toString(0).toString();
	}

	public StringBuilder toString(int tabs) {
		StringBuilder sb=new StringBuilder();
		for (int i=0; i<tabs; i++) {
			sb.append("  ");
		}
		sb.append(decision+" ("+probability+")\n");
		if (left!=null)
			sb.append(left.toString(tabs+1));
		if (right!=null)
			sb.append(right.toString(tabs+1));
		return sb;
	}
}
