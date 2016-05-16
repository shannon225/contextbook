package edu.washington.gs.maccoss.encyclopedia.datastructures;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public enum IonType {
	a,b,c,x,y,z,z1,ap2,bp2,cp2,xp2,yp2,zp2,z1p2,aNL,bNL,cNL,xNL,yNL,zNL,z1NL,ap2NL,bp2NL,cp2NL,xp2NL,yp2NL,zp2NL,z1p2NL;
	
	public static String toString(IonType t) {
		switch (t) {
		case a: return "a";
		case b: return "b";
		case c: return "c";
		case x: return "x";
		case y: return "y";
		case z: return "z";
		case z1: return "z+1";
		case aNL: return "a-NL";
		case bNL: return "b-NL";
		case cNL: return "c-NL";
		case xNL: return "x-NL";
		case yNL: return "y-NL";
		case zNL: return "z-NL";
		case z1NL: return "z+1-NL";

		case ap2: return "a+2H";
		case bp2: return "b+2H";
		case cp2: return "c+2H";
		case xp2: return "x+2H";
		case yp2: return "y+2H";
		case zp2: return "z+2H";
		case z1p2: return "z+1+2H";
		case ap2NL: return "a+2H-NL";
		case bp2NL: return "b+2H-NL";
		case cp2NL: return "c+2H-NL";
		case xp2NL: return "x+2H-NL";
		case yp2NL: return "y+2H-NL";
		case zp2NL: return "z+2H-NL";
		case z1p2NL: return "z+1+2H-NL";
		}
		return "unknown";
	}
	
	@SuppressWarnings("incomplete-switch")
	public static IonType getPlus2(IonType t) {
		switch (t) {
		case a: return ap2;
		case b: return bp2;
		case c: return cp2;
		case x: return xp2;
		case y: return yp2;
		case z: return zp2;
		case z1: return z1p2;
		case aNL: return ap2NL;
		case bNL: return bp2NL;
		case cNL: return cp2NL;
		case xNL: return xp2NL;
		case yNL: return yp2NL;
		case zNL: return zp2NL;
		case z1NL: return z1p2NL;
		}
		throw new EncyclopediaException("Can't make a +2 ion for "+toString(t));
	}
}
