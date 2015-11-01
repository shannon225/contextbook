package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;

public abstract class SQLFile {
	public Connection getConnection(File f) throws IOException {
		Connection c=null;
		try {
			Class.forName("org.sqlite.JDBC");
			c=DriverManager.getConnection("jdbc:sqlite:"+f.getAbsolutePath());
			c.setAutoCommit(false);
			return c;
		} catch (Exception e) {
			System.err.println(e.getClass().getName()+": "+e.getMessage());
			throw new IOException("Error reading database file: "+f.getAbsolutePath());
		}
	}
}
