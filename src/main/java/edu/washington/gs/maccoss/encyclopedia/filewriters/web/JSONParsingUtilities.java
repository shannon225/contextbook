package edu.washington.gs.maccoss.encyclopedia.filewriters.web;

import java.util.ArrayList;
import java.util.List;

public class JSONParsingUtilities {
	public static String removeGroupingBrackets(String value) {
		value = value.trim();
		if (value.startsWith("{") && value.endsWith("}")) {
        	value = value.substring(1, value.length() - 1);
        }
		return value;
	}

	public static String removeArrayBrackets(String value) {
		value = value.trim();
		if (value.startsWith("[") && value.endsWith("]")) {
        	value = value.substring(1, value.length() - 1);
        }
		return value;
	}
    
	public static String getSelectedElement(String[] pairs, String key) {
        for (int i = 0; i < pairs.length; i++) {
        	if (pairs[i].startsWith(key)) {
        		return pairs[i];
        	}
        }
        return null;
    }

	public static String[] splitJsonElements(String json) {
        List<String> elements = new ArrayList<>();
        int braceLevel = 0;
        int bracketLevel = 0;
        int start = 0;

        for (int i = 0; i < json.length(); i++) {
            char ch = json.charAt(i);
            if (ch == '{') braceLevel++;
            else if (ch == '}') braceLevel--;
            else if (ch == '[') bracketLevel++;
            else if (ch == ']') bracketLevel--;
            else if (ch == ',' && braceLevel == 0 && bracketLevel == 0) {
                elements.add(json.substring(start, i));
                start = i + 1;
            }
        }
        elements.add(json.substring(start));

        return elements.toArray(new String[0]);
    }
}
