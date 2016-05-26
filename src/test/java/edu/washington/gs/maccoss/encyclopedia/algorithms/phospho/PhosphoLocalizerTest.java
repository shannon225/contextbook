package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class PhosphoLocalizerTest {
	public static final Color[] colors=new Color[] {Color.red, Color.blue, Color.green, Color.cyan, Color.magenta}; 
	public static void main(String[] args) throws Exception {
		File libraryFile=new File("/Users/searleb/Documents/school/projects/VillenJ_Exactive_HumanPhosphoproteome.elib");
		File diaFile=new File("/Users/searleb/Documents/school/projects/mzml/q06048_rl_MCF7_IMAC_GpX_3.dia");

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		
		LibraryFile library=new LibraryFile();
		library.openFile(libraryFile);
		
		StripeFileInterface stripefile=MzmlToDIAConverter.getFile(diaFile, parameters);
		
		PhosphoLocalizer localizer=new PhosphoLocalizer(stripefile, library, parameters);
		
		double precursorMZ;
		byte precursorCharge;
		String peptideModSeq;
		float retentionTime;
		if (true) {
			precursorMZ=500.730213;
			precursorCharge=(byte)2;
			peptideModSeq="MQS[+80.0]LSLNK";
			retentionTime=1198.3428f;
		} else if (true) {
			precursorMZ=500.899664;
			precursorCharge=(byte)3;
			peptideModSeq="SRPTS[+80.0]FADELAAR";
			retentionTime=1591.183f;
		} else if (true) {
			precursorMZ=503.551374;
			precursorCharge=(byte)3;
			peptideModSeq="A[+42.0]QRHS[+80.0]DSSLEEK";
			retentionTime=517.8737f;
			
		} else {
			precursorMZ=503.272853;
			precursorCharge=(byte)3;
			peptideModSeq="KLS[+80.0]SGDLRVPVTR";
			retentionTime=1309.1414f;
		}
		
		PSMData psmdata=new PSMData(new HashSet<String>(), 0, precursorMZ, precursorCharge, peptideModSeq, retentionTime, 0, 0, 12);
		HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> allVsUniqueList=localizer.runPhosphoLocalization(psmdata, stripefile.getStripes(psmdata.getPrecursorMZ(), 0, Float.MAX_VALUE, false)).getTraces();
		
		ArrayList<Color> shades=new ArrayList<Color>(Arrays.asList(colors));
		ArrayList<XYTrace> traces=new ArrayList<XYTrace>();
		for (Entry<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> entry : allVsUniqueList.entrySet()) {
			String seq=entry.getKey();
			Pair<TFloatFloatHashMap, TFloatFloatHashMap> pair=entry.getValue();
			Color color=shades.remove(0);
			traces.add(new XYTrace(pair.x, GraphType.line, "ALL_"+seq, color, 5.0f));
			traces.add(new XYTrace(pair.y, GraphType.line, "UNI_"+seq, color, 3.0f));
		}
		
		Charter.launchChart("RT", "Score", false, traces.toArray(new XYTrace[traces.size()]));
	}

}
