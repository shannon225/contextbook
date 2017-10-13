package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.File;
import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.math3.stat.inference.TestUtils;

import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ThesaurusElibParser.Coordinate;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ThesaurusElibParser.QuantitationLog;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;

public class EncyclopediaElibParser {
	public static HashMap<String, Coordinate> sampleKey=new HashMap<>();
	private static final int numberOfSampleTypes=6;
	private static final int numberOfReplicates=6;
	private static String getSampleName(int i) {
		switch (i) {
		case 1: return "0 hr";
		case 2: return "2 hr";
		case 3: return "4 hr";
		case 4: return "8 hr";
		case 5: return "16 hr";
		case 6: return "24 hr";
		default: return i > 0 && i < 27 ? String.valueOf((char)(i + 64)) : null;
		}
	}

	public static void loadMap() {
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_1a.mzML", new Coordinate(1, 1));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_1b.mzML", new Coordinate(1, 2));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_1c.mzML", new Coordinate(1, 3));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_1d.mzML", new Coordinate(1, 4));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_1e.mzML", new Coordinate(1, 5));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_1f.mzML", new Coordinate(1, 6));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_2a.mzML", new Coordinate(2, 1));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_2b.mzML", new Coordinate(2, 2));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_2c.mzML", new Coordinate(2, 3));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_2d.mzML", new Coordinate(2, 4));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_2e.mzML", new Coordinate(2, 5));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_2f.mzML", new Coordinate(2, 6));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_3a.mzML", new Coordinate(3, 1));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_3b.mzML", new Coordinate(3, 2));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_3c.mzML", new Coordinate(3, 3));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_3d.mzML", new Coordinate(3, 4));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_3e.mzML", new Coordinate(3, 5));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_3f.mzML", new Coordinate(3, 6));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_4a.mzML", new Coordinate(4, 1));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_4b.mzML", new Coordinate(4, 2));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_4c.mzML", new Coordinate(4, 3));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_4d.mzML", new Coordinate(4, 4));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_4e.mzML", new Coordinate(4, 5));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_4f.mzML", new Coordinate(4, 6));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_5a.mzML", new Coordinate(5, 1));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_5b.mzML", new Coordinate(5, 2));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_5c.mzML", new Coordinate(5, 3));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_5d.mzML", new Coordinate(5, 4));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_5e.mzML", new Coordinate(5, 5));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_5f.mzML", new Coordinate(5, 6));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_6a.mzML", new Coordinate(6, 1));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_6b.mzML", new Coordinate(6, 2));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_6c.mzML", new Coordinate(6, 3));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_6d.mzML", new Coordinate(6, 4));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_6e.mzML", new Coordinate(6, 5));
		sampleKey.put("23aug2017_hela_serum_timecourse_wide_6f.mzML", new Coordinate(6, 6));
	}
	
	/*private static HashMap<String, QuantitationLog> getQuantData(File file) throws IOException, SQLException {
		System.out.println("Parsing "+file.getName()+"...");
		LibraryFile library=new LibraryFile();
		library.openFile(file);

		HashMap<String, QuantitationLog> quantLog=new HashMap<>();
		Connection c=library.getConnection();
		Statement s=c.createStatement();
		ResultSet rs = s.executeQuery("select pep.PrecursorCharge, pep.PeptideModSeq, pep.PeptideSeq, pep.SourceFile, max(pep.LocalizedIntensity), max(pep.TotalIntensity), pep.IsSiteSpecific, pep.RTInSecondsCenter,pep.localizationScore,"+
				"group_concat(p.ProteinAccession, '" + PSMData.ACCESSION_TOKEN + "') as ProteinAccessions " +
				"from " +
				"peptidelocalizations pep " +
				"left join peptidetoprotein p " +
				"where " +
				"pep.PeptideSeq = p.PeptideSeq " +
				"group by pep.rowid;"
		);
		while (rs.next()) {
			Coordinate coord=sampleKey.get(sourceFile);
			if (coord==null) {
				System.out.println("FAILED TO FIND SAMPLE: "+sourceFile);
				System.exit(1);
			}
			QuantitationLog log=quantLog.get(protein);
			if (log==null) {
				log=new QuantitationLog(protein);
				quantLog.put(protein, log);
			}
			log.addIntensity(coord, totalIntensity, rtInSeconds);
		}
		rs.close();
		s.close();
		c.close();
		
		return quantLog;
	}*/
	
	private static double getANOVAPValue(float[][] data) {
		ArrayList<double[]> classes=new ArrayList<>();
		for (int i=0; i<data.length; i++) {
			classes.add(General.toDoubleArray(data[i]));
		}
		return TestUtils.oneWayAnovaPValue(classes);
	}
	
	public static class QuantitationLog {
		final String protein;
		final TObjectFloatHashMap<Coordinate> intensities=new TObjectFloatHashMap<>();
		final TFloatArrayList rtInSecondsList=new TFloatArrayList();
		
		public QuantitationLog(String protein) {
			this.protein=protein;
		}
		
		public int getNumMeasurements() {
			return intensities.size();
		}
		
		public boolean isAtLeastOneCaseFull() {
			float[][] data=getData();
			for (int i=0; i<data.length; i++) {
				boolean full=true;
				for (int j=0; j<data[i].length; j++) {
					if (data[i][j]==0.0f) {
						full=false;
						break;
					}
				}
				if (full) return true;
			}
			return false;
		}
		
		public void addIntensity(Coordinate c, float intensity, float rtInSeconds) {
			intensities.adjustOrPutValue(c, intensity, intensity);
			this.rtInSecondsList.add(rtInSeconds);
		}
		public float[][] getData() {
			float[][] results=new float[numberOfSampleTypes][];
			for (int i=0; i<results.length; i++) {
				results[i]=new float[numberOfReplicates];
			}
			intensities.forEachEntry(new TObjectFloatProcedure<Coordinate>() {
				@Override
				public boolean execute(Coordinate a, float b) {
					results[a.sample-1][a.replicate-1]+=b;
					return true;
				}
			});
			
			return results;
		}
		public float[][] getNormalizedData() {
			float[][] data=getData();
			float[][] normalized=new float[data.length][];
			for (int i=0; i<normalized.length; i++) {
				normalized[i]=new float[data[i].length];
			}

			float grandTotal=0.0f;
			for (int rep=0; rep<data[0].length; rep++) {
				float total=0.0f;
				for (int samp=0; samp<data.length; samp++) {
					total+=data[samp][rep];
				}
				for (int samp=0; samp<data.length; samp++) {
					normalized[samp][rep]=data[samp][rep]/total;
				}
				grandTotal+=total;
			}
			grandTotal=grandTotal/data[0].length;
			for (int rep=0; rep<data[0].length; rep++) {
				for (int samp=0; samp<data.length; samp++) {
					normalized[samp][rep]=normalized[samp][rep]*grandTotal;
				}
			}
			
			return normalized;
		}
	}
	
	public static class Coordinate {
		final int replicate;
		final int sample;
		
		public Coordinate(int replicate, int sample) {
			this.replicate=replicate;
			this.sample=sample;
		}
		
		@Override
		public boolean equals(Object obj) {
			if (obj==null) return false;
			if (!(obj instanceof Coordinate)) return false;
			return hashCode()==obj.hashCode();
		}
		
		@Override
		public int hashCode() {
			return sample+replicate*1000;
		}
	}
}
