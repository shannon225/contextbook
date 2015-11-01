package edu.washington.gs.maccoss.encyclopedia.utils;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.zip.DataFormatException;
import java.util.zip.Deflater;
import java.util.zip.GZIPInputStream;
import java.util.zip.Inflater;

public class CompressionUtils {
	public static byte[] compress(byte[] data) throws IOException {
		Deflater deflater=new Deflater();
		deflater.setInput(data);

		ByteArrayOutputStream outputStream=new ByteArrayOutputStream(data.length);

		deflater.finish();
		byte[] buffer=new byte[1024];
		while (!deflater.finished()) {
			int count=deflater.deflate(buffer);
			outputStream.write(buffer, 0, count);
		}
		outputStream.close();
		byte[] output=outputStream.toByteArray();

		deflater.end();

		return output;
	}

	public static byte[] decompress(byte[] data, int uncompressedLength) throws IOException, DataFormatException {
		Inflater decompresser=new Inflater();
		decompresser.setInput(data);
		byte[] result=new byte[uncompressedLength];
		decompresser.inflate(result);
		decompresser.end();
		return result;
	}

	public static byte[] decompressGzip(byte[] data, int uncompressedLength) throws IOException {
		ByteArrayInputStream bytein=new ByteArrayInputStream(data);
		GZIPInputStream gzin=new GZIPInputStream(bytein);
		ByteArrayOutputStream byteout=new ByteArrayOutputStream();

		int res=0;
		byte[] buf=new byte[1024];
		while (res>=0) {
			res=gzin.read(buf, 0, buf.length);
			if (res>0) {
				byteout.write(buf, 0, res);
			}
		}
		return byteout.toByteArray();
	}
}
