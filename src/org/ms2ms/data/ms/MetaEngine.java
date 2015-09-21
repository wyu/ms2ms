package org.ms2ms.data.ms;

import com.google.common.collect.Multimap;
import com.google.common.collect.Table;
import com.google.common.collect.TreeBasedTable;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import org.ms2ms.alg.LCMSMS;
import org.ms2ms.alg.Peptides;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.util.Map;
import java.util.TreeMap;

/** A concept of meta engine where the scoring attributes are collected from a single or multiple search engines.
 *  The input are one or multiple related MS/MS spectra.
 *
 *
 * Created by Wen Yu on 9/19/2015.
 */
public class MetaEngine
{
  private Table<String, String, SpectrumIdentifier> mRunScanSpectrum    = TreeBasedTable.create();
  private Table<String, Engine, PeptideMatch>       mPeptideEngineMatch = TreeBasedTable.create();

  public MetaEngine() { super(); }

  // Getters and Setters
  public SpectrumIdentifier getEffectiveSpectrumID()
  {
    // TODO need to create a composite spec ID if we are dealing with super consensus of multiple related spectra such as charge variants
    return mRunScanSpectrum!=null?Tools.front(mRunScanSpectrum.values()):null;
  }
  public Table<String, String, SpectrumIdentifier> getSpectrumIDs()        { return mRunScanSpectrum; }
  public Table<String, Engine, PeptideMatch>       getPeptideEngineMatch() { return mPeptideEngineMatch; }

  // add the PSM one at a time
  public MetaEngine add(SpectrumIdentifier ms, PeptideMatch m, Engine engine)
  {
    if (ms!=null && m!=null && engine!=null)
    {
      mRunScanSpectrum.put(Strs.stripLastOf(ms.getSpectrumFile().get(), '.'), LCMSMS.toScanStr(ms.getScanNumbers()), ms);
      String       p = Peptides.toNumModSequence(m);
      PeptideMatch M = mPeptideEngineMatch.get(p, engine);
      // keep only the best one
      if (M==null ||
         (M.hasScore(engine.getCanonicalScore()) &&
          M.getScore(engine.getCanonicalScore())<m.getScore(engine.getCanonicalScore())))
        mPeptideEngineMatch.put(p, engine, m);
    }

    return this;
  }
  // An utility function to align the typical PSMs from a search engines
  public static Map<String, MetaEngine> add(Map<String, MetaEngine> row_matches, Engine engine,
                                            Multimap<SpectrumIdentifier, PeptideMatch> id_match)
  {
    if (row_matches==null) row_matches = new TreeMap<>();

    if (engine!=null && Tools.isSet(id_match))
      for (SpectrumIdentifier id : id_match.keySet())
      {
        String   row = id.getSpectrumFile().get()+"#"+ LCMSMS.toScanStr(id.getScanNumbers());
        MetaEngine E = row_matches.get(row);
        if (E==null) { E = new MetaEngine(); row_matches.put(row, E); }
        // add the contents of the matches
        for (PeptideMatch m : id_match.get(id))
        {
          E.add(id, m, engine);
        }
      }

    return row_matches;
  }
}
