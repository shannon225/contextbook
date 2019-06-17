package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class MaxquantMSMSConverterTest {
	public static void main(String[] args) {
		String[] fileOrder=new String[] {"QE1_24apr2019_BCS_Sample1_DDA_15_0", "QE1_24apr2019_BCS_Sample1_DDA_15_0_2", "QE1_24apr2019_BCS_Sample1_DDA_15_0_3", "QE1_24apr2019_BCS_Sample1_DDA_a_10_5",
				"QE1_24apr2019_BCS_Sample1_DDA_b_7_8", "QE1_24apr2019_BCS_Sample1_DDA_c_4_15", "QE1_24apr2019_BCS_Sample1_DDA_d_3_27", "QE1_24apr2019_BCS_Sample1_DDA_e_2_41",
				"QE1_24apr2019_BCS_Sample1_DDA_f_2_91"};
		
		File tsvFile=new File("/Volumes/searle_ssd/malaria/malaria-dda-maxquant/msms.txt");

		try {
			// protein/file/peptide
			final HashMap<String, HashMap<String, TObjectFloatHashMap<String>>> table=new HashMap<>();
			TableParserMuscle muscle=new TableParserMuscle() {
				@Override
				public void processRow(Map<String, String> row) {
					String rawFile=row.get("Raw file");
					String proteins=row.get("Proteins");
					
					if (!proteins.startsWith("PF")) return; //skip non parasite proteins
					if (proteins.indexOf(";")>0) return; // skip multiprotein peptides

					String precursorKey=row.get("Modified sequence")+"+"+Byte.parseByte(row.get("Charge"));
					float precursorIntensity=Float.parseFloat(row.get("Precursor Intensity"));
					
					HashMap<String, TObjectFloatHashMap<String>> fileMap=table.get(proteins);
					if (fileMap==null) {
						fileMap=new HashMap<>();
						table.put(proteins, fileMap);
					}
					
					TObjectFloatHashMap<String> peptideMap=fileMap.get(rawFile);
					if (peptideMap==null) {
						peptideMap=new TObjectFloatHashMap<>();
						fileMap.put(rawFile, peptideMap);
					}
					
					float previousIntensity=peptideMap.get(precursorKey);
					if (precursorIntensity>previousIntensity) {
						peptideMap.put(precursorKey, precursorIntensity);
					}
				}
				
				@Override
				public void cleanup() {
				}
			};
			
			TableParser.parseTSV(tsvFile, muscle);
			
			System.out.print("Protein");
			for (int i=0; i<fileOrder.length; i++) {
				System.out.print("\t"+fileOrder[i]);
			}
			System.out.println();
			
			for (Entry<String, HashMap<String, TObjectFloatHashMap<String>>> entry : table.entrySet()) {
				String protein=entry.getKey();
				HashMap<String, TObjectFloatHashMap<String>> fileMap=entry.getValue();
				
				System.out.print(protein);
				for (int i=0; i<fileOrder.length; i++) {
					TObjectFloatHashMap<String> peptideMap=fileMap.get(fileOrder[i]);
					float sum=0;
					if (peptideMap!=null) sum=General.sum(peptideMap.values());
					System.out.print("\t"+sum);
				}
				System.out.println();
			}

		} catch (Exception e) {
			Logger.errorLine("Error parsing Maxquant msms.txt:");
			Logger.errorException(e);
			throw new EncyclopediaException(e);
		}
	}

}
