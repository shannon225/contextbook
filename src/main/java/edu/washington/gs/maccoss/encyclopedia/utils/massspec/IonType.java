package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Stroke;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public enum IonType {
	a,b,c,x,y,z,z1,ap2,bp2,cp2,xp2,yp2,zp2,z1p2,aNL,bNL,cNL,xNL,yNL,zNL,z1NL,ap2NL,bp2NL,cp2NL,xp2NL,yp2NL,zp2NL,z1p2NL;
	
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
		}
		return "unknown";
	}
	
	public static String toString(IonType t, byte index) {
		switch (t) {
		case a: return "a"+index;
		case b: return "b"+index;
		case c: return "c"+index;
		case x: return "x"+index;
		case y: return "y"+index;
		case z: return "z"+index;
		case z1: return "z"+index+"+1";
		case aNL: return "a"+index+"-NL";
		case bNL: return "b"+index+"-NL";
		case cNL: return "c"+index+"-NL";
		case xNL: return "x"+index+"-NL";
		case yNL: return "y"+index+"-NL";
		case zNL: return "z"+index+"-NL";
		case z1NL: return "z"+index+"+1-NL";

		case ap2: return "a"+index+"+2H";
		case bp2: return "b"+index+"+2H";
		case cp2: return "c"+index+"+2H";
		case xp2: return "x"+index+"+2H";
		case yp2: return "y"+index+"+2H";
		case zp2: return "z"+index+"+2H";
		case z1p2: return "z"+index+"+1+2H";
		case ap2NL: return "a"+index+"+2H-NL";
		case bp2NL: return "b"+index+"+2H-NL";
		case cp2NL: return "c"+index+"+2H-NL";
		case xp2NL: return "x"+index+"+2H-NL";
		case yp2NL: return "y"+index+"+2H-NL";
		case zp2NL: return "z"+index+"+2H-NL";
		case z1p2NL: return "z"+index+"+1+2H-NL";
		}
		return "unknown";
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
		}
		return missingAnnotationFont;
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
		throw new EncyclopediaException("Can't make a +2 ion for "+toString(t));
	}
}
