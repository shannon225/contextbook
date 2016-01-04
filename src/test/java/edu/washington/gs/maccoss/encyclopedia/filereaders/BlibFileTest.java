package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.IOException;
import java.sql.SQLException;

import junit.framework.TestCase;

public class BlibFileTest extends TestCase {
	public void testCreateTables() throws IOException, SQLException {
		BlibFile file=new BlibFile();
		file.openFile();
		
		assertNotNull(file); // just getting here without errors makes sure the table structure is ok
	}
}
