package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.sql.Connection;

import junit.framework.TestCase;

public class SQLFileTest extends TestCase {
	File blibFile=null;
	BlibFile blib=null;
	@Override
	protected void setUp() throws Exception {
		InputStream is=getClass().getResourceAsStream("/empty.blib");
		blibFile=File.createTempFile("empty", ".blib");
		blibFile.deleteOnExit();
		Files.copy(is, blibFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
		
		blib=new BlibFile();
		blib.openFile(blibFile);
	}

	public void testGetConnection() throws Exception {
		Connection c=blib.getConnection(blibFile);
		assertTrue(c!=null);
		c.close();
	}

	public void testDoesColumnExist() throws Exception {
		assertTrue(blib.doesColumnExist(blibFile, "RefSpectra", "retentionTime"));
	}

}
