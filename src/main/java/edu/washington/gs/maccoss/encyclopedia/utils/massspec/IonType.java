package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Stroke;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;

public enum IonType {
	a,aNL,ap2,ap2NL,
	b,bNL,bp2,bp2NL,
	c,cNL,cp2,cp2NL,
	x,xNL,xp2,xp2NL,
	y,yNL,yp2,yp2NL,
	z,zNL,z1,z1NL,zp2,zp2NL,z1p2,z1p2NL,
	unknown;
	
	private static final String _1="+1";
	private static final String NL="-NL";
	private static final String _1_NL="+1-NL";
	private static final String _2H="+2H";
	private static final String _1_2H="+1+2H";
	private static final String _2H_NL="+2H-NL";
	private static final String _1_2H_NL="+1+2H-NL";
	
	public static final Color oddColor=new Color(26, 148, 49);
	public static final Color bcColor=new Color(226, 75, 59);
	public static final Color yzColor=new Color(59, 109, 226);
	public static final Color missingColor=Color.black;
	public static final BasicStroke primaryStroke=new BasicStroke(2.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
	public static final BasicStroke secondaryStroke=new BasicStroke(1.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
	public static final BasicStroke missingStroke=new BasicStroke(1.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
	public static final Font primaryAnnotationFont=new Font("News Gothic MT", Font.BOLD, 10);
	public static final Font secondaryAnnotationFont=new Font("News Gothic MT", Font.PLAIN, 10);
	public static final Font missingAnnotationFont=new Font("News Gothic MT", Font.PLAIN, 10);

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
		case unknown: return "unknown";
		}
		return "unknown";
	}
	
	public static IonType fromString(String s) {
		if("a".equals(s)) return a;
		if("b".equals(s)) return b;
		if("c".equals(s)) return c;
		if("x".equals(s)) return x;
		if("y".equals(s)) return y;
		if("z".equals(s)) return z;
		if("z+1".equals(s)) return z1;
		if("a-NL".equals(s)) return aNL;
		if("b-NL".equals(s)) return bNL;
		if("c-NL".equals(s)) return cNL;
		if("x-NL".equals(s)) return xNL;
		if("y-NL".equals(s)) return yNL;
		if("z-NL".equals(s)) return zNL;
		if("z+1-NL".equals(s)) return z1NL;

		if("a+2H".equals(s)) return ap2;
		if("b+2H".equals(s)) return bp2;
		if("c+2H".equals(s)) return cp2;
		if("x+2H".equals(s)) return xp2;
		if("y+2H".equals(s)) return yp2;
		if("z+2H".equals(s)) return zp2;
		if("z+1+2H".equals(s)) return z1p2;
		if("a+2H-NL".equals(s)) return ap2NL;
		if("b+2H-NL".equals(s)) return bp2NL;
		if("c+2H-NL".equals(s)) return cp2NL;
		if("x+2H-NL".equals(s)) return xp2NL;
		if("y+2H-NL".equals(s)) return yp2NL;
		if("z+2H-NL".equals(s)) return zp2NL;
		if("z+1+2H-NL".equals(s)) return z1p2NL;

		return unknown;
	}
	
	public static String toString(IonType t, byte index) {
		switch (t) {
		case a: return "a"+index;
		case b: return "b"+index;
		case c: return "c"+index;
		case x: return "x"+index;
		case y: return "y"+index;
		case z: return "z"+index;
		case z1: return "z"+index+_1;
		case aNL: return "a"+index+NL;
		case bNL: return "b"+index+NL;
		case cNL: return "c"+index+NL;
		case xNL: return "x"+index+NL;
		case yNL: return "y"+index+NL;
		case zNL: return "z"+index+NL;
		case z1NL: return "z"+index+_1_NL;

		case ap2: return "a"+index+_2H;
		case bp2: return "b"+index+_2H;
		case cp2: return "c"+index+_2H;
		case xp2: return "x"+index+_2H;
		case yp2: return "y"+index+_2H;
		case zp2: return "z"+index+_2H;
		case z1p2: return "z"+index+_1_2H;
		case ap2NL: return "a"+index+_2H_NL;
		case bp2NL: return "b"+index+_2H_NL;
		case cp2NL: return "c"+index+_2H_NL;
		case xp2NL: return "x"+index+_2H_NL;
		case yp2NL: return "y"+index+_2H_NL;
		case zp2NL: return "z"+index+_2H_NL;
		case z1p2NL: return "z"+index+_1_2H_NL;
		case unknown: return "unknown";
		}
		return "unknown";
	}
	
	public static Pair<IonType, Byte> fromIndexedString(String s) {
		byte type;
		byte index;
		if (s.endsWith(_1_2H_NL)) {
			type=7;
			index=Byte.parseByte(s.substring(1, s.length()-_1_2H_NL.length()));
		} else if (s.endsWith(_2H_NL)) {
			type=6;
			index=Byte.parseByte(s.substring(1, s.length()-_2H_NL.length()));
		} else if (s.endsWith(_1_2H)) {
			type=5;
			index=Byte.parseByte(s.substring(1, s.length()-_1_2H.length()));
		} else if (s.endsWith(_2H)) {
			type=4;
			index=Byte.parseByte(s.substring(1, s.length()-_2H.length()));
		} else if (s.endsWith(_1_NL)) {
			type=3;
			index=Byte.parseByte(s.substring(1, s.length()-_1_NL.length()));
		} else if (s.endsWith(NL)) {
			type=2;
			index=Byte.parseByte(s.substring(1, s.length()-NL.length()));
		} else if (s.endsWith(_1)) {
			type=1;
			index=Byte.parseByte(s.substring(1, s.length()-_1.length()));
		} else {
			type=0;
			index=Byte.parseByte(s.substring(1));
		}
		
		if(s.charAt(0)=='a') {
			switch (type) {
				case 0: return new Pair<IonType, Byte>(a, index); // normal
				case 1: throw new EncyclopediaException("Can't process ion type: "+s); // +1 
				case 2: return new Pair<IonType, Byte>(aNL, index); // neutral loss
				case 3: throw new EncyclopediaException("Can't process ion type: "+s); // +1 neutral loss 
				case 4: return new Pair<IonType, Byte>(ap2, index); // +2H
				case 5: throw new EncyclopediaException("Can't process ion type: "+s); // +1+2H
				case 6: return new Pair<IonType, Byte>(ap2NL, index); // +2H neutral loss
				case 7: throw new EncyclopediaException("Can't process ion type: "+s); // +1+2H neutral loss
			}
		} else if(s.charAt(0)=='b') {
			switch (type) {
				case 0: return new Pair<IonType, Byte>(b, index); // normal
				case 1: throw new EncyclopediaException("Can't process ion type: "+s); // +1 
				case 2: return new Pair<IonType, Byte>(bNL, index); // neutral loss
				case 3: throw new EncyclopediaException("Can't process ion type: "+s); // +1 neutral loss 
				case 4: return new Pair<IonType, Byte>(bp2, index); // +2H
				case 5: throw new EncyclopediaException("Can't process ion type: "+s); // +1+2H
				case 6: return new Pair<IonType, Byte>(bp2NL, index); // +2H neutral loss
				case 7: throw new EncyclopediaException("Can't process ion type: "+s); // +1+2H neutral loss
			}
		} else if(s.charAt(0)=='c') {
			switch (type) {
				case 0: return new Pair<IonType, Byte>(c, index); // normal
				case 1: throw new EncyclopediaException("Can't process ion type: "+s); // +1 
				case 2: return new Pair<IonType, Byte>(cNL, index); // neutral loss
				case 3: throw new EncyclopediaException("Can't process ion type: "+s); // +1 neutral loss 
				case 4: return new Pair<IonType, Byte>(cp2, index); // +2H
				case 5: throw new EncyclopediaException("Can't process ion type: "+s); // +1+2H
				case 6: return new Pair<IonType, Byte>(cp2NL, index); // +2H neutral loss
				case 7: throw new EncyclopediaException("Can't process ion type: "+s); // +1+2H neutral loss
			}
		} else if(s.charAt(0)=='x') {
			switch (type) {
				case 0: return new Pair<IonType, Byte>(x, index); // normal
				case 1: throw new EncyclopediaException("Can't process ion type: "+s); // +1 
				case 2: return new Pair<IonType, Byte>(xNL, index); // neutral loss
				case 3: throw new EncyclopediaException("Can't process ion type: "+s); // +1 neutral loss 
				case 4: return new Pair<IonType, Byte>(xp2, index); // +2H
				case 5: throw new EncyclopediaException("Can't process ion type: "+s); // +1+2H
				case 6: return new Pair<IonType, Byte>(xp2NL, index); // +2H neutral loss
				case 7: throw new EncyclopediaException("Can't process ion type: "+s); // +1+2H neutral loss
			}
		} else if(s.charAt(0)=='y') {
			switch (type) {
				case 0: return new Pair<IonType, Byte>(y, index); // normal
				case 1: throw new EncyclopediaException("Can't process ion type: "+s); // +1 
				case 2: return new Pair<IonType, Byte>(yNL, index); // neutral loss
				case 3: throw new EncyclopediaException("Can't process ion type: "+s); // +1 neutral loss 
				case 4: return new Pair<IonType, Byte>(yp2, index); // +2H
				case 5: throw new EncyclopediaException("Can't process ion type: "+s); // +1+2H
				case 6: return new Pair<IonType, Byte>(yp2NL, index); // +2H neutral loss
				case 7: throw new EncyclopediaException("Can't process ion type: "+s); // +1+2H neutral loss
			}
		} else if(s.charAt(0)=='z') {
			switch (type) {
				case 0: return new Pair<IonType, Byte>(z, index); // normal
				case 1: return new Pair<IonType, Byte>(z1, index); // +1 
				case 2: return new Pair<IonType, Byte>(zNL, index); // neutral loss
				case 3: return new Pair<IonType, Byte>(z1NL, index); // +1 neutral loss 
				case 4: return new Pair<IonType, Byte>(zp2, index); // +2H
				case 5: return new Pair<IonType, Byte>(z1p2, index); // +1+2H
				case 6: return new Pair<IonType, Byte>(zp2NL, index); // +2H neutral loss
				case 7: return new Pair<IonType, Byte>(z1p2NL, index); // +1+2H neutral loss
			}
		}
		
		throw new EncyclopediaException("Can't process ion type: "+s);
	}
	
	public static Color getColor(IonType t) {
		switch (t) {
		case a: return oddColor;
		case b: return bcColor;
		case c: return bcColor;
		case x: return oddColor;
		case y: return yzColor;
		case z: return yzColor;
		case z1: return yzColor;
		case aNL: return oddColor;
		case bNL: return bcColor;
		case cNL: return bcColor;
		case xNL: return oddColor;
		case yNL: return yzColor;
		case zNL: return yzColor;
		case z1NL: return yzColor;

		case ap2: return oddColor;
		case bp2: return bcColor;
		case cp2: return bcColor;
		case xp2: return oddColor;
		case yp2: return yzColor;
		case zp2: return yzColor;
		case z1p2: return yzColor;
		case ap2NL: return oddColor;
		case bp2NL: return bcColor;
		case cp2NL: return bcColor;
		case xp2NL: return oddColor;
		case yp2NL: return yzColor;
		case zp2NL: return yzColor;
		case z1p2NL: return yzColor;
		case unknown: return missingColor;
		}
		return missingColor;
	}
	
	public static Stroke getStroke(IonType t) {
		switch (t) {
		case a: return secondaryStroke;
		case b: return primaryStroke;
		case c: return primaryStroke;
		case x: return secondaryStroke;
		case y: return primaryStroke;
		case z: return secondaryStroke;
		case z1: return primaryStroke;
		case aNL: return secondaryStroke;
		case bNL: return primaryStroke;
		case cNL: return primaryStroke;
		case xNL: return secondaryStroke;
		case yNL: return primaryStroke;
		case zNL: return secondaryStroke;
		case z1NL: return primaryStroke;

		case ap2: return secondaryStroke;
		case bp2: return secondaryStroke;
		case cp2: return secondaryStroke;
		case xp2: return secondaryStroke;
		case yp2: return secondaryStroke;
		case zp2: return secondaryStroke;
		case z1p2: return secondaryStroke;
		case ap2NL: return secondaryStroke;
		case bp2NL: return secondaryStroke;
		case cp2NL: return secondaryStroke;
		case xp2NL: return secondaryStroke;
		case yp2NL: return secondaryStroke;
		case zp2NL: return secondaryStroke;
		case z1p2NL: return secondaryStroke;
		case unknown: return missingStroke;
		}
		return missingStroke;
	}
	
	public static Font getFont(IonType t) {
		switch (t) {
		case a: return secondaryAnnotationFont;
		case b: return primaryAnnotationFont;
		case c: return primaryAnnotationFont;
		case x: return secondaryAnnotationFont;
		case y: return primaryAnnotationFont;
		case z: return secondaryAnnotationFont;
		case z1: return primaryAnnotationFont;
		case aNL: return secondaryAnnotationFont;
		case bNL: return primaryAnnotationFont;
		case cNL: return primaryAnnotationFont;
		case xNL: return secondaryAnnotationFont;
		case yNL: return primaryAnnotationFont;
		case zNL: return secondaryAnnotationFont;
		case z1NL: return primaryAnnotationFont;

		case ap2: return secondaryAnnotationFont;
		case bp2: return secondaryAnnotationFont;
		case cp2: return secondaryAnnotationFont;
		case xp2: return secondaryAnnotationFont;
		case yp2: return secondaryAnnotationFont;
		case zp2: return secondaryAnnotationFont;
		case z1p2: return secondaryAnnotationFont;
		case ap2NL: return secondaryAnnotationFont;
		case bp2NL: return secondaryAnnotationFont;
		case cp2NL: return secondaryAnnotationFont;
		case xp2NL: return secondaryAnnotationFont;
		case yp2NL: return secondaryAnnotationFont;
		case zp2NL: return secondaryAnnotationFont;
		case z1p2NL: return secondaryAnnotationFont;
		case unknown: return missingAnnotationFont;
		}
		return missingAnnotationFont;
	}
	
	public static IonType getCanonicalIonType(IonType t) {
		switch (t) {
		case a: return a;
		case b: return b;
		case c: return c;
		case x: return x;
		case y: return y;
		case z: return z;
		case z1: return z;
		case aNL: return a;
		case bNL: return b;
		case cNL: return c;
		case xNL: return x;
		case yNL: return y;
		case zNL: return z;
		case z1NL: return z;

		case ap2: return a;
		case bp2: return b;
		case cp2: return c;
		case xp2: return x;
		case yp2: return y;
		case zp2: return z;
		case z1p2: return z;
		case ap2NL: return a;
		case bp2NL: return b;
		case cp2NL: return c;
		case xp2NL: return x;
		case yp2NL: return y;
		case zp2NL: return z;
		case z1p2NL: return z;
		case unknown: return unknown;
		}
		return unknown;
	}
	
	@SuppressWarnings("incomplete-switch")
	public static IonType getNL(IonType t) {
		switch (t) {
		case a: return aNL;
		case b: return bNL;
		case c: return cNL;
		case x: return xNL;
		case y: return yNL;
		case z: return zNL;
		case z1: return z1NL;
		case ap2: return ap2NL;
		case bp2: return bp2NL;
		case cp2: return cp2NL;
		case xp2: return xp2NL;
		case yp2: return yp2NL;
		case zp2: return zp2NL;
		case z1p2: return z1p2NL;
		}
		return t; // already a NL
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
		throw new EncyclopediaException("Can't make a +2 ion for type: "+toString(t));
	}
}
