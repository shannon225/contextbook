package edu.washington.gs.maccoss.encyclopedia.utils.io;

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
		try (InputStream is = uri.toURL().openStream()) {
			Files.copy(is, destination, StandardCopyOption.REPLACE_EXISTING);
		}
	}
}
