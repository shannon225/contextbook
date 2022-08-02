package edu.washington.gs.maccoss.encyclopedia.utils;

import java.util.Arrays;

public class StringUtils {
	public static String getCommonName(String[] names, String insert) {
		if (insert==null) insert="";
		if (names==null||names.length<=1) return insert;
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
}
