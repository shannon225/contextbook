package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.RetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.list.array.TFloatArrayList;

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

	
	public static PercolatorData filterData(PercolatorData perc, ArrayList<PeptideScoringResult> data, RetentionTimeFilter filter, SearchParameters parameters) {
		ArrayList<PercolatorPSM> psms=perc.getPsms(); //already reverse sorted

		HashSet<String> passingPSMIDs=new HashSet<String>();
		for (PercolatorPSM psm : psms) {
			if (psm.getQValue()<=parameters.getPercolatorThreshold()) {
				passingPSMIDs.add(psm.getPsmID());
			}
		}
		
		TFloatArrayList deltas=new TFloatArrayList();

		for (PeptideScoringResult result : data) {
			if (result.getGoodStripes().size()>0) {
				String peptideModSeq=result.getEntry().getPeptideModSeq();
				if (passingPSMIDs.contains(peptideModSeq+"+"+result.getEntry().getPrecursorCharge())) {
					LibraryEntry entry=result.getEntry();
					float entryTime=filter.getYValue(entry.getRetentionTime());
					
					Pair<ScoredObject<Stripe>, float[]> first=result.getGoodStripes().get(0);
					float deltaRT=first.x.y.getScanStartTime()/60f-entryTime;
					deltas.add(deltaRT);					
				}
			}
		}
		
		float[] deltaArray=deltas.toArray();
		float mean=General.mean(deltaArray);
		float stdev=General.stdev(deltaArray);
		float upperThreshold=mean+2.0f*stdev;
		float lowerThreshold=mean-2.0f*stdev;

		HashSet<String> rtFilteredPSMIDs=new HashSet<String>();
		for (PeptideScoringResult result : data) {
			if (result.getGoodStripes().size()>0) {
				String peptideModSeq=result.getEntry().getPeptideModSeq();
				LibraryEntry entry=result.getEntry();
				float entryTime=filter.getYValue(entry.getRetentionTime());

				Pair<ScoredObject<Stripe>, float[]> first=result.getGoodStripes().get(0);
				float deltaRT=first.x.y.getScanStartTime()/60f-entryTime;

				if (deltaRT<=upperThreshold&&deltaRT>=lowerThreshold) {
					rtFilteredPSMIDs.add(peptideModSeq);
				}
			}
		}

		int decoys=0;
		int nondecoys=0;
		TFloatArrayList qvalues=new TFloatArrayList();
		for (PercolatorPSM psm : psms) {
			if (rtFilteredPSMIDs.contains(psm.getPsmID())) {
				if (psm.isDecoy()) {
					decoys++;
				} else {
					nondecoys++;
				}
			}
			qvalues.add(decoys/(float)(decoys+nondecoys));
		}
		
		// convert FDRs to q-values
		float minValue=Float.MAX_VALUE;
		for (int i=qvalues.size()-1; i>=0; i--) {
			if (qvalues.get(i)>minValue) {
				qvalues.set(i, minValue);
			} else {
				minValue=qvalues.get(i);
			}
		}

		ArrayList<PercolatorPSM> rtFilteredPSMs=new ArrayList<PercolatorPSM>();
		for (int i=0; i<psms.size(); i++) {
			if (rtFilteredPSMIDs.contains(psms.get(i).getPsmID())) {
				rtFilteredPSMs.add(psms.get(i).clone(qvalues.get(i)));
			}
		}
		
		return perc.clone(rtFilteredPSMs);
	}
}
