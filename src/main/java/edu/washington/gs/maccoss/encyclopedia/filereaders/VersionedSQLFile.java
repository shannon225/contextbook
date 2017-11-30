package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.utils.io.Version;

import java.io.IOException;
import java.sql.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A SQLFile with a metadata table.
 * Version is stored in the metadata table.
 *
 * n.b., All subclasses must call createNewTables() from their open() method
 *
 * Created with IntelliJ IDEA.
 * User: caleb
 * Date: 11/29/2017
 * Time: 11:24 AM
 */
public abstract class VersionedSQLFile extends SQLFile {
	public static final String VERSION_STRING="version";

	/**
	 *
	 * @return The maximal current version of this file. New files will be given this version.
	 */
	public abstract Version getMostRecentVersion();

	protected abstract Connection getConnection() throws IOException, SQLException;

	/**
	 * Apply all alterations to the sql file to update it from its current version to the most recent version
	 * @param currentVersion the current version of the file, or 0.0.0 if it did not previously have a version
	 * @param statement
	 */
	protected abstract void applyPatches(Version currentVersion, Statement statement) throws IOException, SQLException;

	/**
	 * Create the tables as they would appear in the newest version. Any changes
	 * from previous versions should be captured in the patches.
	 * @param statement
	 */
	protected abstract void createTables(Statement statement) throws IOException, SQLException;

	public final void addMetadata(String key, String value) throws IOException, SQLException {
		HashMap<String, String> data=new HashMap<>();
		data.put(key, value);
		addMetadata(data);
	}

	public final void addMetadata(Map<String, String> data) throws IOException, SQLException {
		Connection c=getConnection();
		try {
			PreparedStatement prep=c.prepareStatement("insert or replace into metadata (Key, Value) VALUES (?,?)");
			try {
				for (Map.Entry<String, String> entry : data.entrySet()) {
					prep.setString(1, entry.getKey());
					prep.setString(2, entry.getValue());
					prep.addBatch();
				}
				prep.executeBatch();
				prep.close();
				c.commit();
			} finally {
				prep.close();
			}
		} finally {
			c.close();
		}
	}

	public final void setFileVersion() throws IOException, SQLException {
		addMetadata(VERSION_STRING, getMostRecentVersion().toString());
	}

	public final Version getVersion() throws IOException, SQLException {
		Map<String, String> meta=getMetadata();
		return new Version(meta.get(VERSION_STRING));
	}

	public final HashMap<String, String> getMetadata() throws IOException, SQLException {
		Connection c=getConnection();
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select Key, Value from metadata");

				HashMap<String, String> map=new HashMap<String, String>();
				while (rs.next()) {
					String key=rs.getString(1);
					String value=rs.getString(2);
					map.put(key, value);
				}

				return map;
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

	/**
	 * Create a metadata table if it does not exist. Apply patches to update file, then call the abstract method
	 * to create file-specific tables. This method MUST be called when opening a VersionedSQLFile
	 * @throws IOException
	 * @throws SQLException
	 */
	protected final void createNewTables() throws IOException, SQLException {
		Connection c = getConnection();
		try {
			Statement s=c.createStatement();
			try {
				s.execute("create table if not exists metadata ( Key string not null, Value string not null, primary key (Key) )");
				s.execute("create index if not exists 'Key_Metadata_index' on 'metadata' ('Key' ASC)");

				c.commit(); // commit so getVersion() finds the new table

				Version version = getVersion(); // will be 0.0.0 if metadata table was just created

				applyPatches(version, s);

				createTables(s);

				c.commit();
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}

		setFileVersion();
	}
}
