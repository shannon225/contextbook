package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;

public interface StripeFileInterface {

	HashMap<Range, Float> getRanges();

	void openFile(File userFile) throws IOException, SQLException;

	void openFile() throws IOException, SQLException;

	ArrayList<PrecursorScan> getPrecursors(float minRT, float maxRT) throws IOException, SQLException, DataFormatException;

	ArrayList<Stripe> getStripes(double targetMz, float minRT, float maxRT, boolean sqrt) throws IOException, SQLException;

	ArrayList<Stripe> getStripes(Range targetMzRange, float minRT, float maxRT, final boolean sqrt) throws IOException, SQLException;
	
	float getTIC() throws IOException, SQLException;
	
	float getGradientLength() throws IOException, SQLException;
	
	void close();

	File getFile();
	
	String getOriginalFileName();
}