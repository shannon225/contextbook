package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.net.InetAddress;
import java.net.NetworkInterface;
import java.net.SocketException;
import java.net.UnknownHostException;

public class Networking {
	final protected static char[] hexArray = "0123456789ABCDEF".toCharArray();
	private static final String[] poopAddresses=new String[] {
		"98-90-96-DF-EE-BC","98-5A-EB-D0-2F-5F"
	};
	
	public static void main(String arg[]) {
		try {
			byte[] mac=getMacAddress();
			System.out.println(bytesToHex(mac));
		} catch (Exception e) {
			System.out.println(e);
		}
	}
	
	public static boolean isOffendingAddress() {
		try {
			byte[] mac=getMacAddress();
			String address=bytesToHex(mac);
			for (String b : poopAddresses) {
				if (b.equals(address)) return true;
			}
			return false;
		} catch (Exception e) {
			return false;
		}
	}
	
	public static byte[] getMacAddress() throws UnknownHostException, SocketException {
		InetAddress address=InetAddress.getLocalHost();
		NetworkInterface nwi=NetworkInterface.getByInetAddress(address);
		byte mac[]=nwi.getHardwareAddress();
		return mac;
	}
	
	public static String bytesToHex(byte[] bytes) {
		StringBuilder sb=new StringBuilder();
	    for ( int j = 0; j < bytes.length; j++ ) {
	        int v = bytes[j] & 0xFF;
	        if (sb.length()>0) sb.append('-');
	        sb.append(hexArray[v >>> 4]);
	        sb.append(hexArray[v & 0x0F]);
	    }
	    return sb.toString();
	}
}
