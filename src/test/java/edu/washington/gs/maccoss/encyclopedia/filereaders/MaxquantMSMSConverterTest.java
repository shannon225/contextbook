package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import com.google.common.collect.ImmutableList;
import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.library.EncyclopediaParametersPanel;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.map.hash.TObjectFloatHashMap;
import junit.framework.TestCase;
import org.apache.commons.io.FileUtils;

public class MaxquantMSMSConverterTest extends TestCase {
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

	private static final SearchParameters parameters = new EncyclopediaParametersPanel(new SearchPanel(ProgramType.EncyclopeDIA)).getParameters();

	private Path fasta;
	private Path elib;

	@Override
	protected void setUp() throws Exception {

		fasta = getFasta();
		elib = getEmptyElib();
	}

	@Override
	protected void tearDown() throws Exception {
		FileUtils.deleteQuietly(fasta.toFile());
		fasta = null;

		FileUtils.deleteQuietly(elib.toFile());
		elib = null;
	}

	public void testNewMsmsTxt() throws Exception {
		final Path tsv = getResourceAsFile("msms-new.txt", ".msms.txt");

		final LibraryFile libraryFile = MaxquantMSMSConverter.convertFromMSMSTSV(
				tsv.toFile(),
				fasta.toFile(),
				elib.toFile(),
				parameters
		);

		assertNotNull(libraryFile);
		assertEquals(5, libraryFile.getAllEntries(false, parameters.getAAConstants()).size());
	}

	public void testBadMsmsTxt() throws Exception {
		// A copy of msms-new.txt with a single value from a row's "Intensity" column replaced with an empty string
		final Path tsv = getResourceAsFile("msms-bad.txt", ".msms.txt");

		LibraryFile libraryFile = null;
		try {
			libraryFile = MaxquantMSMSConverter.convertFromMSMSTSV(
					tsv.toFile(),
					fasta.toFile(),
					elib.toFile(),
					parameters
			);
			System.err.println("WARNING! Did not encounter an exception parsing problematic file!");
		} catch (NumberFormatException e) {
			// swallow expected exception and don't run any assertions
			return;
		} catch (Throwable t) {
			throw new AssertionError("Caught unexpected exception!", t);
		}

		assertNotNull(libraryFile);

		// If the parsing doesn't fail, it should include all the entries,
		// or at least all non-erroneous rows (note that the problematic
		// row is neither the first nor the last in the file).
		final int entryCount = libraryFile.getAllEntries(false, parameters.getAAConstants()).size();

		assertTrue("Unexpected number of entries: " + entryCount, ImmutableList.of(4,5).contains(entryCount));
	}

	private static Path getEmptyElib() throws IOException {
		final Path tmp = Files.createTempFile("test_", ".elib");
		tmp.toFile().deleteOnExit();
		return tmp;
	}

	private static Path getFasta() throws IOException {
		return getResourceAsFile("msms-test.fasta", ".fasta");
	}

	private static Path getResourceAsFile(String relativeResourceName, String suffix) throws IOException {
		final Path tsv = Files.createTempFile("test_", suffix);
		try (InputStream tsvResource = MaxquantMSMSConverterTest.class.getResourceAsStream(relativeResourceName)) {
			Files.copy(tsvResource, tsv, StandardCopyOption.REPLACE_EXISTING);
		}
		tsv.toFile().deleteOnExit();
		return tsv;
	}
}
