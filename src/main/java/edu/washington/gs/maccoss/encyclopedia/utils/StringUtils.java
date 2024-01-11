package edu.washington.gs.maccoss.encyclopedia.utils;

import java.util.Arrays;

public class StringUtils {
	public static String getCommonName(String[] names, String insert) {
		if (insert==null) insert="";
		if (names==null||names.length==0) return insert;
		if (names.length==1) return names[0]+insert;
		String standard=names[0];
		
		int pre;
		PREFIX: for (pre=0; pre<standard.length(); pre++) {
			char c=standard.charAt(pre);
			for (int i=1; i<names.length; i++) {
				if (names[i].charAt(pre)!=c) {
					break PREFIX;
				}
			}
		}
		
		if (pre==standard.length()-1) return insert;
		
		int post;
		POSTFIX: for (post=0; post<standard.length(); post++) {
			char c=standard.charAt(standard.length()-1-post);
			for (int i=1; i<names.length; i++) {
				if (names[i].charAt(names[i].length()-1-post)!=c) {
					break POSTFIX;
				}
			}
		}
		
		String preName=names[0].substring(0, pre);
		String postName=names[0].substring(names[0].length()-post);

		return preName+insert+postName;
	}
	public static String[] getUniquePortion(String[] names) {
		if (names==null||names.length<=1) return names;
		String standard=names[0];
		
		int pre;
		PREFIX: for (pre=0; pre<standard.length(); pre++) {
			char c=standard.charAt(pre);
			for (int i=1; i<names.length; i++) {
				if (names[i].charAt(pre)!=c) {
					break PREFIX;
				}
			}
		}
		
		if (pre==standard.length()-1) return names;
		
		int post;
		POSTFIX: for (post=0; post<standard.length(); post++) {
			char c=standard.charAt(standard.length()-1-post);
			for (int i=1; i<names.length; i++) {
				if (names[i].charAt(names[i].length()-1-post)!=c) {
					break POSTFIX;
				}
			}
		}

		String[] newNames=new String[names.length];
		for (int i=0; i<newNames.length; i++) {
			newNames[i]=names[i].substring(pre, names[i].length()-post);
		}
		return newNames;
	}
	
	public static int countSubstring(String s, String sub) {
		int lastIndex = 0;
		int count = 0;
		while (lastIndex != -1) {
			lastIndex = s.indexOf(sub, lastIndex);
			if (lastIndex != -1) {
				count++;
				lastIndex += sub.length();
			}
		}
		return count;
	}
	
	public static int getIndexOf(String[] strings, String target) {
		for (int i=0; i<strings.length; i++) {
			if (target.equals(strings[i])) return i;
		}
		return -1;
	}
	
	public static String getPad(int padLength, char c) {
		char[] ca=new char[padLength];
		Arrays.fill(ca, c);
		return new String(ca);
	}
	
	public static String wrap(String s, int wrap) {
		StringBuilder sb=new StringBuilder();
		char[] ca=s.toCharArray();
		for (int i=0; i<ca.length; i++) {
			if (i>0&&i%wrap==0) sb.append('\n');
			sb.append(ca[i]);
		}
		return sb.toString();
	}
	
	public static String scientificNotation(double v, int precision) {
		assert(precision>=0);
		String s=String.format("%."+precision+"E", v);
		int e=s.indexOf('E');
		String first=s.substring(0, e);
		String second=s.substring(e+1);
		
		int power=Integer.parseInt(second);
		char[] powerChars=Integer.toString(power).toCharArray();
		StringBuilder sb=new StringBuilder(first);
		sb.append("x10");
		for (char c : powerChars) {
			char unicode;
			switch (c) {
			case 0+48: unicode='\u2070'; break;
			case 1+48: unicode='\u00B9'; break;
			case 2+48: unicode='\u00B2'; break;
			case 3+48: unicode='\u00B3'; break;
			case 4+48: unicode='\u2074'; break;
			case 5+48: unicode='\u2075'; break;
			case 6+48: unicode='\u2076'; break;
			case 7+48: unicode='\u2077'; break;
			case 8+48: unicode='\u2078'; break;
			case 9+48: unicode='\u2079'; break;
			default: unicode=c; break;
			}
			sb.append(unicode);
		}
		return sb.toString();
	}
}
