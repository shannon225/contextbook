package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.StringTokenizer;

import org.apache.commons.lang3.text.WordUtils;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filewriters.FastaWriter;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class FastaReaderUtil {
	public static void main3(String[] args) {
		String[] acc= {"ENST00000203616.variant", "ENST00000221347.variant", "ENST00000243300.variant", "ENST00000250617.variant", "ENST00000254765.variant", "ENST00000260359.variant",
				"ENST00000262858.variant", "ENST00000262940.variant", "ENST00000264079.variant", "ENST00000270509.variant", "ENST00000291860.variant", "ENST00000293771.variant",
				"ENST00000301293.variant", "ENST00000304060.variant", "ENST00000304778.variant", "ENST00000305127.variant", "ENST00000310452.variant", "ENST00000313116.variant",
				"ENST00000317103.variant", "ENST00000317543.variant", "ENST00000318504.variant", "ENST00000318529.variant", "ENST00000318991.variant", "ENST00000319340.variant",
				"ENST00000325250.variant", "ENST00000328078.variant", "ENST00000328300.variant", "ENST00000335822.variant", "ENST00000336077.variant", "ENST00000338222.variant",
				"ENST00000342206.variant", "ENST00000342274.variant", "ENST00000342376.variant", "ENST00000342782.variant", "ENST00000350721.variant", "ENST00000355283.variant",
				"ENST00000356413.variant", "ENST00000356606.variant", "ENST00000356790.variant", "ENST00000357242.variant", "ENST00000357544.variant", "ENST00000357991.variant",
				"ENST00000359061.variant", "ENST00000359361.variant", "ENST00000360279.variant", "ENST00000360319.variant", "ENST00000361499.variant", "ENST00000361603.variant",
				"ENST00000369137.variant", "ENST00000369138.variant", "ENST00000369577.variant", "ENST00000369850.variant", "ENST00000370167.variant", "ENST00000370268.variant",
				"ENST00000370270.variant", "ENST00000370401.variant", "ENST00000370479.variant", "ENST00000370620.variant", "ENST00000370622.variant", "ENST00000371372.variant",
				"ENST00000371628.variant", "ENST00000371906.variant", "ENST00000373344.variant", "ENST00000374519.variant", "ENST00000374737.variant", "ENST00000376809.variant",
				"ENST00000376943.variant", "ENST00000377065.variant", "ENST00000377934.variant", "ENST00000378138.variant", "ENST00000378142.variant", "ENST00000378444.variant",
				"ENST00000378455.variant", "ENST00000378657.variant", "ENST00000378982.variant", "ENST00000379806.variant", "ENST00000379869.variant", "ENST00000379873.variant",
				"ENST00000379876.variant", "ENST00000379878.variant", "ENST00000380470.variant", "ENST00000380956.variant", "ENST00000381127.variant", "ENST00000383101.variant",
				"ENST00000395603.variant", "ENST00000396634.variant", "ENST00000396965.variant", "ENST00000397050.variant", "ENST00000397354.variant", "ENST00000398000.variant",
				"ENST00000402510.variant", "ENST00000404681.variant", "ENST00000406854.variant", "ENST00000409324.variant", "ENST00000412866.variant", "ENST00000413795.variant",
				"ENST00000414849.variant", "ENST00000418814.variant", "ENST00000418931.variant", "ENST00000421544.variant", "ENST00000422373.variant", "ENST00000423345.variant",
				"ENST00000424296.variant", "ENST00000426613.variant", "ENST00000432680.variant", "ENST00000442455.variant", "ENST00000446864.variant", "ENST00000449285.variant",
				"ENST00000449970.variant", "ENST00000455586.variant", "ENST00000457062.variant", "ENST00000480796.variant", "ENST00000481285.variant", "ENST00000494664.variant",
				"ENST00000495517.variant", "ENST00000505868.variant", "ENST00000529681.variant", "ENST00000535861.variant", "ENST00000537104.variant", "ENST00000539731.variant",
				"ENST00000559596.variant", "ENST00000560177.variant", "ENST00000560747.variant", "ENST00000586582.variant", "ENST00000586965.variant", "ENST00000590320.variant",
				"ENST00000600128.variant", "ENST00000601739.variant"};
		
		HashSet<String> targets=new HashSet<>();
		for (int i=0; i<acc.length; i++) {
			targets.add(acc[i]);
		}

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		File f=new File("/Volumes/searle_ssd/malaria/hela_specific_database/hela_specific.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);

		for (FastaEntryInterface entry : entries) {
			if (targets.contains(entry.getAccession())) {
				String mutationText=entry.getAnnotation().substring(entry.getAnnotation().lastIndexOf('.')+1);
				StringTokenizer st=new StringTokenizer(mutationText, ";");
				char[] seq=entry.getSequence().toCharArray();
				while (st.hasMoreTokens()) {
					String token=st.nextToken();
					char originalAA=token.charAt(0);
					char replacementAA=token.charAt(token.length()-1);
					int index=Integer.parseInt(token.subSequence(1, token.length()-1).toString());

					//System.out.println(token+" --> "+originalAA+"|"+index+"|"+replacementAA+" ("+seq[index-1]+")");
					seq[index-1]=originalAA;
				}
				
				String newAnn=entry.getAnnotation().substring(0, entry.getAnnotation().lastIndexOf(' '));
				newAnn=newAnn.replaceAll(".variant", "");
				System.out.print('>');
				System.out.println(newAnn);
				System.out.println(WordUtils.wrap(new String(seq), 80, "\n", true));
				System.out.println();
			}
		}
	}
	public static void main(String[] args) throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();

		Range mzRange=new Range(396.43, 1002.70);
		AminoAcidConstants constants=parameters.getAAConstants();
		File f=new File("/Volumes/searle_ssd/malaria/hela_specific_database/hela_specific.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);
		
		HashSet<String> allPeptides=new HashSet<>();
		DigestionEnzyme enzyme=DigestionEnzyme.getEnzyme("trypsin");
		for (FastaEntryInterface entry : entries) {
			if (entry.getAccession().indexOf(".variant")>=0) {
				continue;
			}
			ArrayList<FastaPeptideEntry> peptides=enzyme.digestProtein(entry, 7, 30, 1, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
			for (FastaPeptideEntry pep : peptides) {
				for (int pepCharge : new int[] {2, 3}) {
					String seq=pep.getSequence();
					double pepMass=constants.getMass(seq)+MassConstants.oh2;
					double pepChargedMass=(pepMass+MassConstants.protonMass*pepCharge)/pepCharge;

					if (mzRange.contains(pepChargedMass)) {
						if (seq.indexOf('B')>=0||seq.indexOf('J')>=0||seq.indexOf('O')>=0||seq.indexOf('U')>=0||seq.indexOf('X')>=0||seq.indexOf('Z')>=0||seq.indexOf('*')>=0) {
							continue;
						} else {
							allPeptides.add(seq);
						}
					}
				}
			}
		}
		HashMap<String, String> variantPeptides=new HashMap<>();
		for (FastaEntryInterface entry : entries) {
			if (entry.getAccession().indexOf(".variant")==-1) {
				continue;
			}
			ArrayList<FastaPeptideEntry> peptides=enzyme.digestProtein(entry, 7, 30, 1, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()), false);
			for (FastaPeptideEntry pep : peptides) {
				for (int pepCharge : new int[] {2, 3}) {
					String seq=pep.getSequence();
					double pepMass=constants.getMass(seq)+MassConstants.oh2;
					double pepChargedMass=(pepMass+MassConstants.protonMass*pepCharge)/pepCharge;

					if (mzRange.contains(pepChargedMass)) {
						if (seq.indexOf('B')>=0||seq.indexOf('J')>=0||seq.indexOf('O')>=0||seq.indexOf('U')>=0||seq.indexOf('X')>=0||seq.indexOf('Z')>=0||seq.indexOf('*')>=0) {
							continue;
						} else {
							if (!allPeptides.contains(seq)) {
								//System.out.println(seq+"\t"+entry.getAnnotation());
								variantPeptides.put(seq, entry.getAnnotation());
							}
							allPeptides.add(seq);
						}
					}
				}
			}
		}
		
		LibraryFile library=new LibraryFile();
		//library.openFile(new File("/Volumes/searle_ssd/malaria/hela_specific_database/narrow/whole_hela_specific_clib.z2_nce33.elib"));
		//library.openFile(new File("/Volumes/searle_ssd/malaria/hela_specific_database/wide_clib/23aug2017_hela_serum_timecourse_pool_wide_001_170829031834.dia.elib"));
		//library.openFile(new File("/Volumes/searle_ssd/malaria/hela_specific_database/wide_clib/23aug2017_hela_serum_timecourse_pool_wide_002.dia.elib"));
		library.openFile(new File("/Volumes/searle_ssd/malaria/hela_specific_database/wide_clib/23aug2017_hela_serum_timecourse_pool_wide_003.dia.elib"));
		
		int count=0;
		ArrayList<LibraryEntry> found=library.getAllEntries(false, new AminoAcidConstants());
		for (LibraryEntry entry : found) {
			if (variantPeptides.containsKey(entry.getPeptideSeq())) {
				System.out.println(entry.getPeptideSeq()+"\t"+entry.getPrecursorMZ()+"\t"+entry.getRetentionTime()+"\t"+variantPeptides.get(entry.getPeptideSeq()));
				count++;
			}
		}
		library.close();
		System.out.println("found "+count);
	}
	public static void main2(String[] args) throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		File f=new File("/Volumes/searle_ssd/malaria/hela_specific_database/HeLa_Database.txt.fasta");
		File fout=new File("/Volumes/searle_ssd/malaria/hela_specific_database/hela_specific.fasta");
		
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(f, parameters);
		
		FastaWriter writer=new FastaWriter(fout);
		
		for (FastaEntryInterface entry : entries) {
			if (entry.getAnnotation().startsWith("mz")) {
				String[] values=entry.getAnnotation().split("\\|");
				if (values.length==3) {
					String accession=values[1];
					String description=values[2];
					String[] descriptionValues=description.split(" ");
					if (descriptionValues[descriptionValues.length-1].startsWith("MU=")) {
						accession=accession+".variant";
					}
					String newAnnotation=accession+" "+description;
					System.out.println(newAnnotation);
					FastaEntry outEntry=new FastaEntry(fout.getName(), newAnnotation, entry.getSequence());
					writer.write(outEntry);
				}
			}
		}
		
		writer.close();
	}
}
