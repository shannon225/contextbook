package edu.washington.gs.maccoss.encyclopedia.utils.io;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorVersion;
import org.junit.Assume;
import org.junit.Ignore;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.net.URI;
import java.net.URLConnection;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import static org.junit.Assert.*;
import static org.junit.Assume.assumeTrue;

public class UriDownloaderTest {
	@Rule
	public TemporaryFolder tmpDir = new TemporaryFolder();

	@Test
	public void downloadFromLocalUri() throws Exception {
		// Get a local file -- which one doesn't really matter.
		final Path source = PercolatorVersion.v3p05.getPercolator().toPath();
		assumeTrue(Files.exists(source));

		final Path destination = tmpDir.newFile().toPath();
		assumeTrue(Files.deleteIfExists(destination));

		final URI uri = source.toUri();

		UriDownloader.downloadFromUri(uri, destination);

		assertTrue(Files.exists(destination));
		assertEquals(Files.size(source), Files.size(destination));
	}

	@Test
	public void downloadFromHttpUri() throws Exception {
		final URI uri = new URI("https://bitbucket.org/searleb/encyclopedia/downloads/encyclopedia_logo.png");

		final URLConnection conn = uri.toURL().openConnection();
		conn.connect();

		final long size = conn.getContentLength();

		final Path destination = tmpDir.newFile().toPath();
		Files.deleteIfExists(destination);

		UriDownloader.downloadFromUri(uri, destination);

		assertTrue(Files.exists(destination));
		assertEquals(size, Files.size(destination));
	}

	@Test
	@Ignore
	public void downloadFromFtpUri() throws Exception {
		final URI uri = new URI(""); //TODO: test FTP URI

		final URLConnection conn = uri.toURL().openConnection();
		conn.connect();

		final long size = conn.getContentLength();

		final Path destination = tmpDir.newFile().toPath();
		Files.deleteIfExists(destination);

		UriDownloader.downloadFromUri(uri, destination);

		assertTrue(Files.exists(destination));
		assertEquals(size, Files.size(destination));
	}
}