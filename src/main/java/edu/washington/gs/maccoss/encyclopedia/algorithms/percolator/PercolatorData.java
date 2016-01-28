package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class PercolatorData {
	private final int majorVersion;
	private final int minorVersion;
	private final String percolatorVersion;
	private final String commandLine;
	private final String otherCommandLine;
	private final float pi0PSMs;
	private final float pi0Peptides;
	private final int psmsQLevel;
	private final int peptidesQLevel;
	private final ArrayList<PercolatorPSM> psms;
	
	public PercolatorData(int majorVersion, int minorVersion, String percolatorVersion, String commandLine, String otherCommandLine, float pi0PSMs, float pi0Peptides, int psmsQLevel,
			int peptidesQLevel, ArrayList<PercolatorPSM> psms) {
		this.majorVersion=majorVersion;
		this.minorVersion=minorVersion;
		this.percolatorVersion=percolatorVersion;
		this.commandLine=commandLine;
		this.otherCommandLine=otherCommandLine;
		this.pi0PSMs=pi0PSMs;
		this.pi0Peptides=pi0Peptides;
		this.psmsQLevel=psmsQLevel;
		this.peptidesQLevel=peptidesQLevel;
		this.psms=new ArrayList<PercolatorPSM>(psms);
		Collections.sort(this.psms);
		Collections.reverse(this.psms);
	}
	
	
	
	public PercolatorData clone(ArrayList<PercolatorPSM> newPSMs) {
		return new PercolatorData(majorVersion, minorVersion, percolatorVersion, commandLine, otherCommandLine, pi0PSMs, pi0Peptides, psmsQLevel, peptidesQLevel, newPSMs);
	}



	public ArrayList<PercolatorPSM> getPsms() {
		return psms;
	}
	
	public void writeToFile(File outputFile) {
		PrintWriter writer=null;

		try {
			writer=new PrintWriter(outputFile, "UTF-8");

			writer.println("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"+
					"<percolator_output \n"+
					"xmlns=\"http://per-colator.com/percolator_out/14\" \n"+
					"xmlns:p=\"http://per-colator.com/percolator_out/14\" \n"+
					"xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" \n"+
					"xsi:schemaLocation=\"http://per-colator.com/percolator_out/14 https://github.com/percolator/percolator/raw/pout-1-4/src/xml/percolator_out.xsd\" \n"+
					"p:majorVersion=\""+majorVersion+"\" p:minorVersion=\""+minorVersion+"\" p:percolator_version=\""+percolatorVersion+"\">\n");
			
			writer.println(" <process_info>\n"+
					"    <command_line>"+commandLine+"</command_line>\n"+
					"    <other_command_line>"+otherCommandLine+"</other_command_line>\n"+
					"    <pi_0_psms>"+pi0PSMs+"</pi_0_psms>\n"+
					"    <pi_0_peptides>"+pi0Peptides+"</pi_0_peptides>\n"+
					"    <psms_qlevel>"+psmsQLevel+"</psms_qlevel>\n"+
					"    <peptides_qlevel>"+peptidesQLevel+"</peptides_qlevel>\n"+
					"  </process_info>\n");
			
			writer.println("  <psms>");
			for (PercolatorPSM percolatorPSM : psms) {
				writer.print(percolatorPSM.toString());
			}
			writer.println("  </psms>");
			
			// FIXME write peptides
			writer.flush();
		} catch (FileNotFoundException e) {
			throw new EncyclopediaException("Error setting up output file: "+outputFile.getAbsolutePath(), e);
		} catch (UnsupportedEncodingException e) {
			throw new EncyclopediaException("Error setting up output file: "+outputFile.getAbsolutePath(), e);
		} finally {
			writer.close();
		}
	}
}
