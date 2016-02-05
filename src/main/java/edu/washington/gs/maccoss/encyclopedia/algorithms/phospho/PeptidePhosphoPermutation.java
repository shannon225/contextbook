package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

public class PeptidePhosphoPermutation {
	private static final String phospho = "[+79.96633]";
	public static final String[] charges = new String[] { "+", "++", "+++",
			"++++", "+5", "+6" };

	public static void main(String[] args) {
		HashSet<String> alreadySeen = new HashSet<String>();

		for (int i = 0; i < PeptideList.peptideList.length; i++) {
			String peptide = PeptideList.peptideList[i];
			if (alreadySeen.contains(peptide)) {
				continue;
			} else {
				alreadySeen.add(peptide);
			}

			int charge = PeptideList.chargeList[i];
			double rt=PeptideList.rtList[i];

			int count = 0;
			StringBuilder sb = new StringBuilder();
			CharacterIterator iterator = new CharacterIterator(peptide);
			ArrayList<Integer> stys = new ArrayList<Integer>();
			while (iterator.hasNext()) {
				char c = iterator.next();

				if ('[' == c) {
					count++;
					while (true) {
						c = iterator.next();
						if (']' == c) {
							break;
						}
					}
				} else {
					if ('S' == c || 'T' == c || 'Y' == c) {
						stys.add(sb.length());
					}
					sb.append(c);
				}
			}

			if (count == 1) {
				int previous=0;
				for (int j = 0; j < stys.size(); j++) {
					StringBuilder seq = new StringBuilder(sb);
					Integer thisIndex = stys.get(j);
					seq.insert(thisIndex + 1, phospho);
					
					String annotatedSeq=new String(seq).replaceAll("C", "C[+57.021464]");
					//System.out.print(annotatedSeq+"\t");
					System.out.print(annotatedSeq+charges[charge-1]);

					System.out.print("\t"+rt);
					
					System.out.print("\tb");
					for (int m = previous; m < thisIndex; m++) {
						if (m>previous) System.out.print(',');
						System.out.print((m+1));
					}
					System.out.print("\ty");
					for (int m = previous; m < thisIndex; m++) {
						if (m>previous) System.out.print(',');
						System.out.print((sb.length()-m-1));
					}
					System.out.println();
					
					previous=thisIndex;
				}
			} else if (count ==2) {
				int previous=0;
				for (int j = 0; j < stys.size(); j++) {
					Integer thisIndex = stys.get(j); // FIXME only prints out localization for the first site!
					for (int k = j+1; k < stys.size(); k++) {
						StringBuilder seq = new StringBuilder(sb);
						seq.insert(stys.get(j) + 1, phospho);
						seq.insert(stys.get(k) + 1 + phospho.length(), phospho);
						
						String annotatedSeq=new String(seq).replaceAll("C", "C[+57.021464]");
						//System.out.print(annotatedSeq+"\t");
						System.out.print(annotatedSeq+charges[charge-1]);

						System.out.print("\t"+rt);
						
						System.out.print("\tb");
						for (int m = previous; m < thisIndex; m++) {
							if (m>previous) System.out.print(',');
							System.out.print((m+1));
						}
						System.out.print("\ty");
						for (int m = previous; m < thisIndex; m++) {
							if (m>previous) System.out.print(',');
							System.out.print((sb.length()-m-1));
						}
						System.out.println();
						
						previous=thisIndex;
					}
				}
			} else {
				throw new RuntimeException("No support for >2 phosphos!");
			}
		}
	}

	static class CharacterIterator implements Iterator<Character> {

		private final String str;
		private int pos = 0;

		public CharacterIterator(String str) {
			this.str = str;
		}

		public boolean hasNext() {
			return pos < str.length();
		}

		public Character next() {
			return str.charAt(pos++);
		}

		public void remove() {
			throw new UnsupportedOperationException();
		}
	}
}
