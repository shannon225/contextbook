package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.IOException;
import java.lang.ref.WeakReference;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;

/**
 * This cache keeps a fixed number of entire stripes in memory. If a request is
 * made for an unloaded stripe, then it is added to the hardCache. After enough
 * values enter the hardCache, it is truncated and the remainder are added to
 * the softCache.
 * 
 * CACHE IS NOT THREAD SAFE! ONLY USE TO SET UP THREADED JOBS!
 *
 */
public class WeakReferenceStripeCache {
	private final int hardCacheSize;
	
	// hardCache is bounded to hardCacheSize
	private final LinkedList<Pair<Range, ArrayList<FragmentScan>>> hardCache=new LinkedList<Pair<Range, ArrayList<FragmentScan>>>();
	
	// softCache is unbounded and only cleaned up by promotion requests. This means that there needs to be a relatively low number of Ranges in total
	private final LinkedList<Pair<Range, WeakReference<ArrayList<FragmentScan>>>> softCache=new LinkedList<Pair<Range, WeakReference<ArrayList<FragmentScan>>>>();

	private final ArrayList<Range> allRanges=new ArrayList<Range>();
	protected final StripeFileInterface stripeFile;
	protected final SearchParameters parameters;

	public WeakReferenceStripeCache(StripeFileInterface stripeFile, SearchParameters parameters) {
		this(6, stripeFile, parameters);
	}
	
	public WeakReferenceStripeCache(int hardCacheSize, StripeFileInterface stripeFile, SearchParameters parameters) {
		this.hardCacheSize=hardCacheSize;
		this.stripeFile=stripeFile;
		this.parameters=parameters;
		for (Range range : stripeFile.getRanges().keySet()) {
			if (!parameters.useTargetWindowCenter()||range.contains(parameters.getTargetWindowCenter())) {
				allRanges.add(range);
			}
		}
		Collections.sort(allRanges);
	}
	
	public ArrayList<FragmentScan> getStripes(float mz) {
		// first check hard references
		Pair<Range, ArrayList<FragmentScan>> foundHardPair=null;
		for (Pair<Range, ArrayList<FragmentScan>> pair : hardCache) {
			if (pair.x.contains(mz)) {
				foundHardPair=pair;
				break;
			}
		}
		// if found stripes as hard references, return them!
		if (foundHardPair!=null) {
			hardCache.remove(foundHardPair);
			hardCache.addFirst(foundHardPair);
			return foundHardPair.y;
		}
		
		// then check soft references
		Pair<Range, WeakReference<ArrayList<FragmentScan>>> foundSoftPair=null;
		ArrayList<FragmentScan> foundPreviouslySoftStripes=null; // can be null even if foundSoftPair isn't!
		for (Pair<Range, WeakReference<ArrayList<FragmentScan>>> pair : softCache) {
			if (pair.x.contains(mz)) {
				foundPreviouslySoftStripes=pair.y.get();
				foundSoftPair=pair;
				break;
			}
		}

		Range foundPreviouslySoftRange=null;
		if (foundSoftPair!=null) {
			// if found soft reference, then remove it from the soft queue!
			foundPreviouslySoftRange=foundSoftPair.x;
			softCache.remove(foundSoftPair);
		} else {
			// if not found, then find the original range
			for (Range range : allRanges) {
				if (range.contains(mz)) {
					foundPreviouslySoftRange=range;
					break;
				}
			}
			if (foundPreviouslySoftRange==null) {
				// request outside of intended range!
				return new ArrayList<FragmentScan>();
			}
		}
		
		// if soft reference is null (or wasn't in the cache in the first place), then retrieve
		if (foundPreviouslySoftStripes==null) {
			foundPreviouslySoftStripes=getStripesFromFile(mz);
		}
		
		// insert into hard cache and return
		hardCache.addFirst(new Pair<Range, ArrayList<FragmentScan>>(foundPreviouslySoftRange, foundPreviouslySoftStripes));
		while (hardCache.size()>hardCacheSize) {
			Pair<Range, ArrayList<FragmentScan>> previouslyHardPair=hardCache.removeLast();
			softCache.addFirst(new Pair<Range, WeakReference<ArrayList<FragmentScan>>>(previouslyHardPair.x, new WeakReference<ArrayList<FragmentScan>>(previouslyHardPair.y)));
		}
		return foundPreviouslySoftStripes;
	}
	
	protected ArrayList<FragmentScan> getStripesFromFile(float mz) {
		try {
			ArrayList<FragmentScan> stripes=stripeFile.getStripes(mz, -Float.MAX_VALUE, Float.MAX_VALUE, false);
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
