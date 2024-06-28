package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;

import java.util.ArrayList;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

public enum EncodedAminoAcid {
	A(0),
	C(1),
	D(2),
	E(3),
	F(4),
	G(5),
	H(6),
	I(7),
	K(8),
	L(9),
	M(10),
	N(11),
	P(12),
	Q(13),
	R(14),
	S(15),
	T(16),
	V(17),
	W(18),
	Y(19),
	CCam(20),
	Mox(21),
	Wox(22),
	PyroCCam(23),
	PyroGlu(24),
	SPhos(25),
	TPhos(26),
	YPhos(27),
	KAc(28),
	KSucc(29),
	KUb(30),
	KMe(31),
	KDiMe(32),
	KTriMe(33),
	RMe(34),
	RDiMe(35),
	KTMT0(36),
	KTMT10(37),
	nAc(38),
	nTMT0(39),
	nTMT10(40);
	
	public static int MAX_ENCODING_LENGTH=45; 
	
	private final int index;
	private EncodedAminoAcid(int index) {
		this.index=index;
	}
	
	public int getIndex() {
		return index;
	}
	
	public String toString() {
		switch (this) {
			case A: return "A";
			case C: return "C";
			case D: return "D";
			case E: return "E";
			case F: return "F";
			case G: return "G";
			case H: return "H";
			case I: return "I";
			case K: return "K";
			case L: return "L";
			case M: return "M";
			case N: return "N";
			case P: return "P";
			case Q: return "Q";
			case R: return "R";
			case S: return "S";
			case T: return "T";
			case V: return "V";
			case W: return "W";
			case Y: return "Y";
			case CCam: return "C[+57]";
			case Mox: return "M[+16]";
			case Wox: return "W[+16]";
			case PyroCCam: return "C[+40]";
			case PyroGlu: return "Q[-17]";
			case SPhos: return "S[+80]";
			case TPhos: return "T[+80]";
			case YPhos: return "Y[+80]";
			case KAc: return "K[+42]";
			case KSucc: return "K[+100]";
			case KUb: return "K[+114]";
			case KMe: return "K[+14]";
			case KDiMe: return "K[+28]";
			case KTriMe: return "K[+42]";
			case RMe: return "R[+14]";
			case RDiMe: return "K[+28]";
			case KTMT0: return "K[+224]";
			case KTMT10: return "K[+229]";
			case nAc: return "[+42]";
			case nTMT0: return "[+224]";
			case nTMT10: return "[+229]";
			default: throw new EncyclopediaException("Unexpected amino acid ["+getIndex()+"]!");
		}
	}
	
	public boolean isNTerm() {
		switch (this) {
			case nAc: return true;
			case nTMT0: return true;
			case nTMT10: return true;
			default: return false;
		}
	}

	/**
	 * NOTE: will return NULL if more AAs than maxPeptideLength
	 * @param sequence
	 * @param aminoAcidConstants
	 * @param maxPeptideLength
	 * @return
	 */
    public static INDArray encode(String sequence, AminoAcidConstants aminoAcidConstants, int maxPeptideLength) {
    	EncodedAminoAcid[] aas=EncodedAminoAcid.getAAs(sequence, aminoAcidConstants);
    	if (aas.length>maxPeptideLength) return null;
    	
        INDArray encoded = Nd4j.zeros(maxPeptideLength, MAX_ENCODING_LENGTH);
        
        int start=aas[0].isNTerm()?0:1;
        for (int i = start; i < aas.length; i++) {
            encoded.putScalar(new int[]{i, aas[i].index}, 1.0);
        }
        return encoded.reshape(1, maxPeptideLength * MAX_ENCODING_LENGTH);
    }
	
	public static EncodedAminoAcid[] getAAs(String sequence, AminoAcidConstants aminoAcidConstants) {
		char[] ca=sequence.toCharArray();
		
		ArrayList<EncodedAminoAcid> aas=new ArrayList<EncodedAminoAcid>();
		int i=0;
		while (i<ca.length) {
			char c='n';
			double mass=0.0;
			
			if (isUppercaseAlphabetic(ca[i])) {
				c=ca[i];
				i++;
			}
			
			if (i<ca.length&&ca[i]=='[') {
				StringBuilder sb=new StringBuilder();
				i++;
				while (ca[i]!=']') {
					sb.append(ca[i]);
					i++;
				}
				mass=Double.valueOf(sb.toString());
				i++;
			}
			if (mass!=0.0) {
				mass=aminoAcidConstants.getAccurateModificationMass(c, mass);
			}
			
			if (c=='n') {
				aas.add(getNTermMod(mass));
			
			} else {
				aas.add(getAA(c, mass));
			}
		}
		return aas.toArray(new EncodedAminoAcid[0]);
	}

    private static boolean isUppercaseAlphabetic(char c) {
        return c >= 'A' && c <= 'Z';
    }

	private static final MassTolerance tolerance=new MassTolerance(1.0); // 1 ppm is about the accuracy of floats
	public static EncodedAminoAcid getNTermMod(double mass) {
		if (tolerance.equals(mass, 42.010565)) return nAc;
		if (tolerance.equals(mass, 224.152478)) return nTMT0;
		if (tolerance.equals(mass, 229.162932)) return nTMT10;
		
		throw new EncyclopediaException("Unexpected n-term modification ["+mass+"]!");
	}
	
	public static EncodedAminoAcid getAA(char aa, double mass) {
		if (mass==0.0) {
			switch (aa) {
				case 'A': return A;
				case 'C': return C;
				case 'D': return D;
				case 'E': return E;
				case 'F': return F;
				case 'G': return G;
				case 'H': return H;
				case 'I': return I;
				case 'K': return K;
				case 'L': return L;
				case 'M': return M;
				case 'N': return N;
				case 'P': return P;
				case 'Q': return Q;
				case 'R': return R;
				case 'S': return S;
				case 'T': return T;
				case 'V': return V;
				case 'W': return W;
				case 'Y': return Y;
				
				default: throw new EncyclopediaException("Unexpected amino acid "+aa+"["+mass+"]!");
			}
			
		} else {
			// try to find amino acid, otherwise fall back on encoding for unmodified AA (supports unknown silac options)

			switch (aa) {
				case 'C': 
					if (tolerance.equals(mass, 57.0214635)) return CCam;
					if (tolerance.equals(mass, 57.0214635 - 17.026549)) return PyroCCam;
					return C;
					
				case 'E':
					if (tolerance.equals(mass, -18.010565)) return PyroGlu;
					return E;
				
				case 'K':
					if (tolerance.equals(mass, 42.010565)) return KAc; 
					if (tolerance.equals(mass, 100.016044)) return KSucc; 
					if (tolerance.equals(mass, 114.042927)) return KUb; 
					if (tolerance.equals(mass, 14.015650)) return KMe; 
					if (tolerance.equals(mass, 28.031300)) return KDiMe; 
					if (tolerance.equals(mass, 42.046950)) return KTriMe; 
					if (tolerance.equals(mass, 224.152478)) return KTMT0; 
					if (tolerance.equals(mass, 229.162932)) return KTMT10; 
					return K;
					
				case 'M': 
					if (tolerance.equals(mass, 15.994915)) return Mox;
					return M;
					
				case 'Q': 
					if (tolerance.equals(mass, -17.026549)) return PyroGlu;
					return Q;
					
				case 'R': 
					if (tolerance.equals(mass, 14.015650)) return RMe; 
					if (tolerance.equals(mass, 28.031300)) return RDiMe; 
					return R;
					
				case 'S': 
					if (tolerance.equals(mass, 79.966331)) return SPhos;
					return S;
					
				case 'T': 
					if (tolerance.equals(mass, 79.966331)) return TPhos;
					return T;
					
				case 'W': 
					if (tolerance.equals(mass, 15.994915)) return Wox;
					return W;
					
				case 'Y': 
					if (tolerance.equals(mass, 79.966331)) return YPhos;
					return Y;
					
				case 'F': return F;
				case 'G': return G;
				case 'H': return H;
				case 'I': return I;
				case 'L': return L;
				case 'N': return N;
				case 'P': return P;
				case 'V': return V;
					
				default: throw new EncyclopediaException("Unexpected amino acid "+aa+"["+mass+"]!");
			}
		}
	}
}
