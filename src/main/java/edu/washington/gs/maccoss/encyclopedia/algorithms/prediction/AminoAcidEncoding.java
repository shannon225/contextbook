package edu.washington.gs.maccoss.encyclopedia.algorithms.prediction;

import java.util.ArrayList;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

public enum AminoAcidEncoding {
	// can't start at 0, that's for empty values!
	A(1),
	C(2),
	D(3),
	E(4),
	F(5),
	G(6),
	H(7),
	I(8),
	K(9),
	L(10),
	M(11),
	N(12),
	P(13),
	Q(14),
	R(15),
	S(16),
	T(17),
	V(18),
	W(19),
	Y(20),
	CCam(21),
	Mox(22),
	Wox(23),
	PyroCCam(24),
	PyroGlu(25),
	SPhos(26),
	TPhos(27),
	YPhos(28),
	KAc(29),
	KSucc(30),
	KUb(31),
	KMe(32),
	KDiMe(33),
	KTriMe(34),
	RMe(35),
	RDiMe(36),
	KTMT0(37),
	KTMT10(38),
	nAc(39),
	nTMT0(40),
	nTMT10(41),
	n(42),
	c(43);
	
	public static int MAX_ENCODING_LENGTH=45; 
	
	private final int index;
	private AminoAcidEncoding(int index) {
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
			case n: return "";
			case c: return "";
			default: throw new EncyclopediaException("Unexpected amino acid ["+getIndex()+"]!");
		}
	}
	
	public boolean isNTerm() {
		switch (this) {
			case nAc: return true;
			case nTMT0: return true;
			case nTMT10: return true;
			case n: return true;
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
    	AminoAcidEncoding[] aas=AminoAcidEncoding.getAAs(sequence, aminoAcidConstants);
    	if (aas.length>maxPeptideLength) return null;
    	
        INDArray encoded = Nd4j.zeros(maxPeptideLength, MAX_ENCODING_LENGTH);
        
        int start=aas[0].isNTerm()?0:1;
        for (int i = start; i < aas.length; i++) {
            encoded.putScalar(new int[]{i, aas[i].index}, 1.0);
        }
        return encoded.reshape(1, maxPeptideLength * MAX_ENCODING_LENGTH);
    }
	
	public static AminoAcidEncoding[] getAAs(String sequence, AminoAcidConstants aminoAcidConstants) {
		char[] ca=sequence.toCharArray();
		
		ArrayList<AminoAcidEncoding> aas=new ArrayList<AminoAcidEncoding>();
		aas.add(n);
		
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
				aas.clear(); // remove previous 'n'
				aas.add(getNTermMod(mass));
			
			} else {
				aas.add(getAA(c, mass));
			}
		}
		aas.add(c);
		return aas.toArray(new AminoAcidEncoding[0]);
	}

    private static boolean isUppercaseAlphabetic(char c) {
        return c >= 'A' && c <= 'Z';
    }

	private static final MassTolerance tolerance=new MassTolerance(1.0); // 1 ppm is about the accuracy of floats
	public static AminoAcidEncoding getNTermMod(double mass) {
		if (tolerance.equals(mass, 42.010565)) return nAc;
		if (tolerance.equals(mass, 224.152478)) return nTMT0;
		if (tolerance.equals(mass, 229.162932)) return nTMT10;
		
		throw new EncyclopediaException("Unexpected n-term modification ["+mass+"]!");
	}
	
	public static AminoAcidEncoding getAA(char aa, double mass) {
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
				case 'n': return n;
				case 'c': return c;
				
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
