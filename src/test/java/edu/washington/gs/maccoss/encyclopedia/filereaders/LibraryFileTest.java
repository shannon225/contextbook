package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;

public class LibraryFileTest {
	public static void main(String[] args) throws Exception {
		File f=new File("/Volumes/BriansSSD/ForBrian/CSF_Concensus_Isoform_AD_CN_Library.elib");
		LibraryFile library=(LibraryFile)BlibToLibraryConverter.getFile(f);
		Connection c=library.getConnection();
		try {
			Statement s=c.createStatement();
			try {
				ResultSet rs=s.executeQuery("select PeptideSeq, ProteinAccessions from proteins");

				while (rs.next()) {
					String peptide=rs.getString(1);
					String accessions=rs.getString(2);

					HashSet<String> accessionSet=PSMData.stringToAccessions(accessions);
					ArrayList<String> accessionList=new ArrayList<String>(accessionSet);
					Collections.sort(accessionList);
					System.out.println(peptide+"\t"+PSMData.accessionsToString(accessionList));
				}

			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}
}
