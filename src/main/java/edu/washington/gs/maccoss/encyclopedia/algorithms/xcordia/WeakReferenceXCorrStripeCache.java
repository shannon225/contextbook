package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.WeakReferenceStripeCache;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

/**
 * This cache keeps a fixed number of entire stripes in memory. If a request is
 * made for an unloaded stripe, then it is added to the hardCache. After enough
 * values enter the hardCache, it is truncated and the remainder are added to
 * the softCache.
 * 
 * CACHE IS NOT THREAD SAFE! ONLY USE TO SET UP THREADED JOBS!
 *
 */
public class WeakReferenceXCorrStripeCache extends WeakReferenceStripeCache {
	
	public WeakReferenceXCorrStripeCache(int hardCacheSize, StripeFileInterface stripeFile, SearchParameters parameters) {
		super(hardCacheSize, stripeFile, parameters);
	}

	public WeakReferenceXCorrStripeCache(StripeFileInterface stripeFile, SearchParameters parameters) {
		super(stripeFile, parameters);
	}

	protected ArrayList<FragmentScan> getStripesFromFile(float mz) {
		try {
			ArrayList<FragmentScan> stripes=stripeFile.getStripes(mz, -Float.MAX_VALUE, Float.MAX_VALUE, false);
			ArrayList<FragmentScan> xcorrStripes=new ArrayList<FragmentScan>();
			for (FragmentScan stripe : stripes) {
				xcorrStripes.add(new XCorrStripe(stripe, parameters));
			}
			Collections.sort(stripes);
			return stripes;
			
		} catch (IOException ioe) {
			Logger.errorLine("Error processing "+stripeFile.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", ioe);
		} catch (SQLException sqle) {
			Logger.errorLine("Error processing "+stripeFile.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", sqle);
		}
	}
}
