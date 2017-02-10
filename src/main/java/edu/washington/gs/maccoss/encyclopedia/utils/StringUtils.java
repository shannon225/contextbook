package edu.washington.gs.maccoss.encyclopedia.utils;

public class StringUtils {
	public static String[] getUniquePortion(String[] names) {
		if (names==null||names.length<=1) return names;
		String standard=names[0];
		
		int pre;
		PREFIX: for (pre=0; pre<standard.length(); pre++) {
			char c=standard.charAt(pre);
			for (int i=0; i<names.length; i++) {
				if (names[i].charAt(pre)!=c) {
					break PREFIX;
				}
			}
		}
		
		if (pre==standard.length()-1) return names;
		int post;
		POSTFIX: for (post=standard.length()-1; post>=0; post--) {
			char c=standard.charAt(post);
			for (int i=0; i<names.length; i++) {
				if (names[i].charAt(post)!=c) {
					break POSTFIX;
				}
			}
		}

		String[] newNames=new String[names.length];
		for (int i=0; i<newNames.length; i++) {
			newNames[i]=names[i].substring(pre, post+1);
		}
		return newNames;
	}
}
