package edu.washington.gs.maccoss.encyclopedia.utils;

public class OSDetector {
	final private static String sys=System.getProperty("os.name").toLowerCase();
	final private static String arch=System.getProperty("os.arch").toLowerCase();

	public static void main(String[] args) {
		System.out.println(sys);
		System.out.println(arch);

		System.out.println(getOSName(getOS()));
	}

	public static OS getOS() {
		if (isWindows()) {
			return OS.WINDOWS;
		} else if (isMac()) {
			return OS.MAC;
		} else if (isUnix()) {
			return OS.LINUX;
		}
		throw new EncyclopediaException("Unknown operating system ["+sys+"]");
	}
	
	public static boolean isARM() {
		if ("aarch64".equals(arch)||"arm".equals(arch)) {
			return true;
		}
		return false;
	}
	
	public static String getOSName(OS os) {
		switch (os) {
			case WINDOWS: return "Windows";
			case MAC: return "MacOS";
			case LINUX: return "Unix/Linux";
		}
		return "Unknown";
	}

	public static boolean isWindows() {
		return sys.indexOf("win")>=0;
	}

	public static boolean isMac() {
		return sys.indexOf("mac")>=0;
	}

	public static boolean isUnix() {
		return (sys.indexOf("nix")>=0||sys.indexOf("nux")>=0||sys.indexOf("aix")>0);
	}

	public enum OS {
		WINDOWS,MAC,LINUX;
	}
}
