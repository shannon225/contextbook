package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.zip.DataFormatException;

import gnu.trove.map.hash.TObjectFloatHashMap;

public class IRTdbFile extends SQLFile {
	private final File tempFile;

	public IRTdbFile(File file) {
		this.tempFile=file;
	}

	
	public TObjectFloatHashMap<String> getIRTs() throws IOException, SQLException, DataFormatException {
		TObjectFloatHashMap<String> map=new TObjectFloatHashMap<String>();
		
		Connection c=getConnection(tempFile);
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select PeptideModSeq, Irt from IrtLibrary");
				while (rs.next()) {
					String peptideModSeq=rs.getString(1);
					float retentionTime=(float)rs.getDouble(2);
					map.put(peptideModSeq, retentionTime);
				}
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
		
		return map;
	}

}
