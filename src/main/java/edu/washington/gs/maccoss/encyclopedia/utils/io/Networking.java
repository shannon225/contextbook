package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.net.InetAddress;
import java.net.NetworkInterface;
import java.net.SocketException;
import java.net.UnknownHostException;

public class Networking {
	final protected static char[] hexArray = "0123456789ABCDEF".toCharArray();
	private static final String[] poopAddresses=new String[] {
		"98-90-96-DF-EE-BC",
		"B8-CA-3A-98-6D-BF",
		"F8-32-E4-BD-D1-46"
	};
	
	public static void main(String arg[]) {
		try {
			byte[] mac=getMacAddress();
			System.out.println(bytesToHex(mac));
		} catch (Exception e) {
			System.out.println(e);
		}
	}
	
	public static String getUserID() {
		try {
			byte[] mac=getMacAddress();
			return bytesToHex(mac);

		} catch (Exception e) {
			return "unknown";
		}
	}
	
	public static int isOffendingAddress() {
		try {
			byte[] mac=getMacAddress();
			String address=bytesToHex(mac);
			for (int i=0; i<poopAddresses.length; i++) {
				if (poopAddresses[i].equals(address)) return (i+1);
			}
			return 0;
		} catch (Exception e) {
			return 0;
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
