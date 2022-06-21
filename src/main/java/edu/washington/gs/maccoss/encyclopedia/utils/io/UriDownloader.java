package edu.washington.gs.maccoss.encyclopedia.utils.io;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import org.apache.commons.lang3.NotImplementedException;

import java.io.IOException;
import java.io.InputStream;
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;

/**
 * Utilities for fetching files from a location specified as a URI.
 */
public class UriDownloader {
	/**
	 * Download the contents of {@code uri} to a local file located at {@code destination}.
	 *
	 * Supports any of the following URI schemes:
	 * <ul>
	 *     <li>HTTP/S</li>
	 *     <li>FTP</li>
	 *     <li>Local file ({@code file://})</li>
	 *     <li>AWS S3 ({@code s3://})</li> -- TODO
	 * </ul>
	 */
	public static void downloadFromUri(URI uri, Path destination) throws IOException {
		if ("s3".equalsIgnoreCase(uri.getScheme())) {
			final Pair<String, String> s3coords = parseS3Uri(uri);
			downloadFromS3(s3coords.x, s3coords.y, destination);
			return;
		}

		try (InputStream is = uri.toURL().openStream()) {
			Files.copy(is, destination, StandardCopyOption.REPLACE_EXISTING);
		}
	}

	/**
	 * Download the specified S3 object to the given destination.
	 */
	public static void downloadFromS3(String bucket, String key, Path destination) {
		throw new NotImplementedException("TODO"); //TODO
	}

	/**
	 * Parse the given URI into S3 object coordinates (bucket and key).
	 *
	 * @param uri A URI with the {@code s3://} scheme
	 *
	 * @return A {@link Pair} of the S3 bucket and object key.
	 */
	static Pair<String, String> parseS3Uri(URI uri) {
		if (!"s3".equals(uri.getScheme())) {
			throw new IllegalArgumentException("parseS3Uri() requires a URI with the s3:// scheme!");
		}

		final String bucket, key;
		bucket = uri.getHost();
		key = uri.getPath().substring(1); // Strip leading slash

		return new Pair<>(bucket, key);
	}
}
