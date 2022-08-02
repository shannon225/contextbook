package edu.washington.gs.maccoss.encyclopedia.utils.io;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorVersion;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import org.junit.*;
import org.junit.rules.TemporaryFolder;
import software.amazon.awssdk.core.exception.SdkClientException;
import software.amazon.awssdk.services.s3.S3Client;
import software.amazon.awssdk.services.s3.model.HeadObjectRequest;
import software.amazon.awssdk.services.s3.model.HeadObjectResponse;

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

	@Test
	@Ignore
	public void downloadFromS3Uri() throws Exception {
		final URI uri = new URI("s3://bucket/key"); //TODO: test S3 URI

		final Pair<String, String> pair = UriDownloader.parseS3Uri(uri);
		final String bucket = pair.x, key = pair.y;

		final S3Client s3Client = S3Client.builder().build();
		final HeadObjectRequest headObjectRequest =
				HeadObjectRequest.builder()
						.bucket(bucket)
						.key(key)
						.build();

		final HeadObjectResponse headObjectResponse;
		try {
			headObjectResponse = s3Client.headObject(headObjectRequest);
		} catch (SdkClientException e) {
			// Just ignore this test because we can't access the URI.
			// Most common cause will be lack of credentials, but if we
			// can't access the URI (e.g. wrong account or doesn't exist)
			// we should also ignore the test.

			// Log the issue for clarity
			Logger.errorLine("Skipping S3 download test: " + e.getMessage());

			throw new AssumptionViolatedException("Skipping S3 download test", e);
		}

		final long size = headObjectResponse.contentLength();

		final Path destination = tmpDir.newFile().toPath();
		Files.deleteIfExists(destination);

		UriDownloader.downloadFromUri(uri, destination);

		assertTrue(Files.exists(destination));
		assertEquals(size, Files.size(destination));
	}
}