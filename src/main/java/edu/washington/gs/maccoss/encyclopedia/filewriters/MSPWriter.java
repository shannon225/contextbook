package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PTMMap;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PTMMap.PostTranslationalModification;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;

public class MSPWriter {
	public static final DecimalFormat df=new DecimalFormat("0.0000");
    
	public static void writeMSP(File mspFile, LibraryFile library, SearchParameters parameters) {
		PrintWriter writer=null;
		try {
			writer=new PrintWriter(mspFile);
		for (LibraryEntry rawEntry : library.getAllEntries(false, parameters.getAAConstants())) {
			AnnotatedLibraryEntry entry=new AnnotatedLibraryEntry(rawEntry, parameters);
			
			writer.println("Name: "+entry.getPeptideSeq()+"/"+entry.getPrecursorCharge());
			writer.println("MW: "+MassConstants.getPeptideMass(entry.getPrecursorMZ(), entry.getPrecursorCharge()));
			writer.println("PrecursorMZ: "+entry.getPrecursorMZ());
			writer.println("Charge: "+entry.getPrecursorCharge());
			writer.println("RetentionTimeMins: "+entry.getScanStartTime()/60f);
			
			StringBuilder comment=new StringBuilder("Comment:");
			comment.append(" Parent="+entry.getPrecursorMZ());
			
			FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
			String[] aas=model.getAas();
			double[] deltaMasses=model.getModificationMasses();
			ArrayList<String> ptmStrings=new ArrayList<>(); 
			for (int i = 0; i < deltaMasses.length; i++) {
				if (deltaMasses[i]!=0.0) {
					char c=FragmentationModel.parseAA(aas[i]).x;
					PostTranslationalModification ptm=PTMMap.getPTM(deltaMasses[i], Character.toString(c));
					
					// note, indicies start at 0
					ptmStrings.add(i+","+c+","+ptm.getName());
				}
			}
			
			comment.append(" Mods="+ptmStrings.size()); // reports 0 if empty, annotates using "/" if not
			for (String ptm : ptmStrings) {
				comment.append("/"+ptm);
			}
			
			writer.println(comment.toString());
			
			writer.println("Num peaks: "+entry.getIonCount());
			FragmentIon[] ions=entry.getIonAnnotations();
			float[] intensities=entry.getIntensityArray();
			double[] masses=entry.getMassArray();
			for (int i = 0; i < masses.length; i++) {
				if (intensities[i]>LibraryEntry.minimumIntensityThreshold) {
					String annotation=toAnnotationString(ions[i]);
					if (annotation==null) {
						annotation="?";
					} else {
						double delta=ions[i].getMass()-masses[i];
						
						annotation=annotation+"/"+df.format(delta);
					}
					writer.println(masses[i]+"\t"+intensities[i]+"\t"+annotation);
				}
			}
			
			writer.println();
		}
		
		} catch (FileNotFoundException e) {
			throw new EncyclopediaException("Error parsing library", e);
		} catch (IOException e) {
			throw new EncyclopediaException("Error parsing library", e);
		} catch (SQLException e) {
			throw new EncyclopediaException("Error parsing library", e);
		} catch (DataFormatException e) {
			throw new EncyclopediaException("Error parsing library", e);
		} finally {
			if (writer!=null) {
				writer.close();
			}
		}
	}

	private static String _2H="^2";
	private static String _3H="^3";
	private static String _4H="^4";
	@SuppressWarnings("incomplete-switch")
	public static String toAnnotationString(FragmentIon ion) {
		if (ion==null) return null;
		
		int index=ion.getIndex();
		switch (ion.getType()) {
		case a: return "a"+index;
		case b: return "b"+index;
		case c: return "c"+index;
		case x: return "x"+index;
		case y: return "y"+index;
		case z: return "z"+index;
		case z1: return "z"+index+"+1";

		case ap2: return "a"+index+_2H;
		case bp2: return "b"+index+_2H;
		case cp2: return "c"+index+_2H;
		case xp2: return "x"+index+_2H;
		case yp2: return "y"+index+_2H;
		case zp2: return "z"+index+_2H;
		
		case bp3: return "b"+index+_3H;
		case bp4: return "b"+index+_4H;
		case yp3: return "y"+index+_3H;
		case yp4: return "y"+index+_4H;
		}
		
		return null;
	}
}
