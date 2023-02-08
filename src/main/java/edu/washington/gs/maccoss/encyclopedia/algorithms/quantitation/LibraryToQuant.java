package edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation;

import java.io.File;
import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.precursor.DDAPrecursorIntegrator;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.SubProgressIndicator;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;

public class LibraryToQuant {

	public static void saveQuantData(LibraryFile library, File resultFile, File rawDirectory, boolean integratePrecursors, SearchParameters params, ProgressIndicator progress) {
		try {
			HashMap<String, File> sourceToRawFile = getRawFileMapping(rawDirectory, library);
			saveQuantData(library, resultFile, sourceToRawFile, integratePrecursors, params, progress);
			
		} catch (DataFormatException ioe) {
			Logger.errorLine("Error processing "+library.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", ioe);
		} catch (IOException ioe) {
			Logger.errorLine("Error processing "+library.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", ioe);
		} catch (SQLException ioe) {
			Logger.errorLine("Error processing "+library.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", ioe);
		}
	}
	public static void saveQuantData(LibraryFile library, File resultFile, HashMap<String, File> sourceToRawFile, boolean integratePrecursors, SearchParameters params, ProgressIndicator progress) {
		Pair<TObjectFloatHashMap<String>, ArrayList<IntegratedLibraryEntry>> pair=extractQuantData(library, sourceToRawFile, integratePrecursors, params, progress);

		try {
			LibraryFile elib=new LibraryFile();
			elib.openFile();
			elib.dropIndices();
			
			pair.x.forEachEntry(new TObjectFloatProcedure<String>() {
				@Override
				public boolean execute(String a, float b) {
					try {
						elib.addTIC(a, b);
					} catch (Exception e) {
						throw new EncyclopediaException(e);
					}
					return true;
				}
			});
			elib.addIntegratedEntries(!integratePrecursors, pair.y, Optional.empty(), Optional.empty(), params.getAAConstants(), params.getPercolatorThreshold());

			elib.createIndices();
			elib.saveAsFile(resultFile);
			elib.close();

		} catch (IOException ioe) {
			Logger.errorLine("Error processing "+library.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", ioe);
		} catch (SQLException ioe) {
			Logger.errorLine("Error processing "+library.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", ioe);
		}
		
	}
	public static Pair<TObjectFloatHashMap<String>, ArrayList<IntegratedLibraryEntry>> extractQuantData(LibraryFile library, File rawDirectory, boolean integratePrecursors, SearchParameters params, ProgressIndicator progress) {
		try {
			HashMap<String, File> sourceToRawFile = getRawFileMapping(rawDirectory, library);
			return extractQuantData(library, sourceToRawFile, integratePrecursors, params, progress);
			
		} catch (DataFormatException ioe) {
			Logger.errorLine("Error processing "+library.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", ioe);
		} catch (IOException ioe) {
			Logger.errorLine("Error processing "+library.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", ioe);
		} catch (SQLException ioe) {
			Logger.errorLine("Error processing "+library.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", ioe);
		}
	}
	public static Pair<TObjectFloatHashMap<String>, ArrayList<IntegratedLibraryEntry>> extractQuantData(LibraryFile library, HashMap<String, File> sourceToRawFile, boolean integratePrecursors, SearchParameters params, ProgressIndicator progress) {
		try {
			TObjectFloatHashMap<String> ticMap=new TObjectFloatHashMap<>();
			ArrayList<LibraryEntry> entries=library.getAllEntries(false, params.getAAConstants());
			HashMap<String, ArrayList<LibraryEntry>> sourceSpecificEntries=new HashMap<>();
			for (LibraryEntry entry : entries) {
				String source=entry.getSource();
				ArrayList<LibraryEntry> thisList=sourceSpecificEntries.get(source);
				if (thisList==null) {
					thisList=new ArrayList<>();
					sourceSpecificEntries.put(source, thisList);
				}
				thisList.add(entry);
			}

			ArrayList<IntegratedLibraryEntry> allExtractPeptides=new ArrayList<>();
			for (Map.Entry<String, ArrayList<LibraryEntry>> mapEntry : sourceSpecificEntries.entrySet()) {
				String key = mapEntry.getKey();
				ArrayList<LibraryEntry> entryList = mapEntry.getValue();
				
				ArrayList<PSMData> data=new ArrayList<PSMData>();
				for (LibraryEntry entry : entryList) {
					PSMData psmdata=new PSMData(entry.getAccessions(), entry.getSpectrumIndex(), entry.getPrecursorMZ(), entry.getPrecursorCharge(), 
							entry.getPeptideModSeq(), entry.getRetentionTime(), entry.getScore(), 1.0f-entry.getScore(), params.getExpectedPeakWidth(), false, 
							params.getAAConstants());
					data.add(psmdata);
				}
				
				File rawFile=sourceToRawFile.get(key);
				if (rawFile==null) continue;
				
				StripeFileInterface diaFile=StripeFileGenerator.getFile(rawFile, params);
				ticMap.put(diaFile.getOriginalFileName(), diaFile.getTIC());

				ArrayList<IntegratedLibraryEntry> extractPeptides;
				SubProgressIndicator subProgress=new SubProgressIndicator(progress, 1.0f/sourceToRawFile.size());
				if (integratePrecursors) {
					extractPeptides=DDAPrecursorIntegrator.integrateSearch(subProgress, data, diaFile, params);
				} else {
					PeptideQuantExtractor extractor=new PeptideQuantExtractor(subProgress, diaFile, params);
					extractPeptides=extractor.extractPeptides(data, Optional.ofNullable(null), true);
				}

				Logger.logLine("Attempted extraction for: "+entryList.size()+", found "+extractPeptides.size());
				
				if (entryList.size()!=extractPeptides.size()) {
					HashSet<String> peptideModSeqs=new HashSet<String>();
					for (LibraryEntry entry : extractPeptides) {
						peptideModSeqs.add(entry.getPeptideModSeq());
					}
					for (LibraryEntry entry : entryList) {
						if (!peptideModSeqs.contains(entry.getPeptideModSeq())) {
							Logger.logLine("Failed to extract "+entry.getPeptideModSeq());
						}
					}
				}
				allExtractPeptides.addAll(extractPeptides);
			}
			return new Pair<>(ticMap, allExtractPeptides);

		} catch (InterruptedException ie) {
			Logger.errorLine("Interruption processing "+library.getFile().getName());
			throw new EncyclopediaException("Interruption parsing Stripe file", ie);
		} catch (DataFormatException ioe) {
			Logger.errorLine("Error processing "+library.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", ioe);
		} catch (IOException ioe) {
			Logger.errorLine("Error processing "+library.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", ioe);
		} catch (SQLException ioe) {
			Logger.errorLine("Error processing "+library.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", ioe);
		}
	}
	
	private static HashMap<String, File> getRawFileMapping(File rawDirectory, LibraryFile library) throws IOException, SQLException, DataFormatException {
		Connection c=library.getConnection();
		try {
			Statement s=c.createStatement();
			try {
				ArrayList<String> sourceFiles=new ArrayList<String>();
				
				Logger.logLine("Getting source files...");
				ResultSet rs=s.executeQuery("select distinct SourceFile from entries");
				while (rs.next()) {
					sourceFiles.add(rs.getString(1));
				}
				rs.close();
				
				File[] possibleFiles=rawDirectory.listFiles();
				HashMap<String, File> sourceToRawFile=new HashMap<>();
				
				// warn about missing source files
				for (String source : sourceFiles) {
					File match=null;
					for (File file : possibleFiles) {
						if (file.getName().equalsIgnoreCase(source)) {
							match=file;
							break;
						}
					}
					if (match==null) {
						Logger.errorLine("Cannot find file associated with ["+source+"] in ["+rawDirectory.getAbsolutePath()+"]. Skipping integrating these peptides!");
					} else {
						sourceToRawFile.put(source, match);
					}
				}
				return sourceToRawFile;
			} finally {
				s.close();
			}
		} finally {
			c.close();
		}
	}

}
