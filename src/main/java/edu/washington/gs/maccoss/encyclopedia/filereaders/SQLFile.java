package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public abstract class SQLFile {
	public Connection getConnection(File f) throws IOException {
		Connection c=null;
		try {
			Class.forName("org.sqlite.JDBC");
			c=DriverManager.getConnection("jdbc:sqlite:"+f.getAbsolutePath());
			c.setAutoCommit(false);
			return c;
		} catch (Exception e) {
			Logger.errorLine(e.getClass().getName()+": "+e.getMessage());
			throw new IOException("Error reading database file: "+f.getAbsolutePath());
		}
	}
	
	public boolean doesColumnExist(File f, String table, String column) throws IOException, SQLException {
		Connection c=getConnection(f);
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("SELECT sql FROM sqlite_master WHERE type = 'table' AND name = '"+table+"'");
				while (rs.next()) {
					String statement=rs.getString(1);
					if (statement.toLowerCase().indexOf(" "+column.toLowerCase()+" ")>=0) {
						return true;
					}
				}
				return false;
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}
}
