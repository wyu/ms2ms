package org.ms2ms.data.ms;

import com.google.common.collect.*;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import org.ms2ms.algo.LCMSMS;
import org.ms2ms.algo.PSMs;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.util.*;

/** A concept of meta engine where the scoring attributes are collected from a single or multiple search engines.
 *  The input are one or multiple related MS/MS spectra.
 *
 *
 * Created by Wen Yu on 9/19/2015.
 */
public class MetaPeptideMatch
{
  private Table<String, String, SpectrumIdentifier> mRunScanSpectrum    = TreeBasedTable.create();
  private Table<String, Engine, PeptideMatch>       mPeptideEngineMatch = TreeBasedTable.create();
  private TreeBasedTable<Integer, Double, String>   mVoteRankSeqIndex   = TreeBasedTable.create();

  public MetaPeptideMatch() { super(); }

  // get the best match out of the consensus
  public PeptideMatch getTopPeptideMatch(Map<Engine, PeptideMatch> dirt)
  {
    PeptideMatch top = null;
    if (Tools.isSet(mPeptideEngineMatch))
    {
      if (!Tools.isSet(mVoteRankSeqIndex)) rank();
      // grab the top vote
      Integer   vote = mVoteRankSeqIndex.rowKeySet().first();
      Double    rank = mVoteRankSeqIndex.row(vote).firstKey();
      String peptide = mVoteRankSeqIndex.get(vote, rank);
      // make the match
      top = PSMs.clone(Tools.front(mPeptideEngineMatch.row(peptide).values()));
      // merge the scores
      top.getScoreMap().clear();
      for (Engine E : mPeptideEngineMatch.columnKeySet())
      {
        // fill in a dummy match if not presence already
        if (mPeptideEngineMatch.get(peptide, E)==null)
        {
//          PeptideMatch dirt = null;
//          for (PeptideMatch d : mPeptideEngineMatch.column(E).values())
//            if (dirt==null || dirt.getRank()<d.getRank()) dirt=d;
//          // look for the lowest rank for this engine
//          PeptideMatch dummy = PSMs.clone(dirt); dummy.getScoreMap().clear(); dummy.setRejected(true);
//          // half the scores
//          for (String score : dirt.getScoreMap().keySet())
//            dummy.addScore(score, dirt.getScore(score)/2d);
          mPeptideEngineMatch.put(peptide, E, dirt.get(E));
        }
      }
      double voted=0d;
      for (Engine E : mPeptideEngineMatch.columnKeySet())
      {
        PeptideMatch m = mPeptideEngineMatch.get(peptide, E);
        if (!m.isRejected()) voted+=1d;
        for (String scr : m.getScoreMap().keySet())
        {
          top.addScore(scr, m.getScore(scr));
        }
      }
      top.addScore("Vote", voted); top.setRank(1);
    }

    return top;
  }
  // get the consensus matches or the best single match if no consensus can be found
  public Collection<PeptideMatch> getTopPeptideMatches(Map<Engine, PeptideMatch> dirt)
  {
    if (!Tools.isSet(mPeptideEngineMatch)) return null;

    Collection<PeptideMatch>  tops = new ArrayList<>();

//    if (Tools.front(mPeptideEngineMatch.values()).getNeutralPeptideMass()>1000.4 &&
//        Tools.front(mPeptideEngineMatch.values()).getNeutralPeptideMass()<1000.5)
//    {
////      System.out.println();
//    }
    // looking for the consensus where 2 or more engines voted for the same peptide
    for (String peptide : mPeptideEngineMatch.rowKeySet())
      if (mPeptideEngineMatch.row(peptide).keySet().size()>1)
        tops.add(getCombinedMatch(peptide, dirt));

    // if no consensus was found, pick the top ranked from each engine
//    if (tops.size()==0)
//    {
      for (Engine E : mPeptideEngineMatch.columnKeySet())
        for (PeptideMatch m : mPeptideEngineMatch.column(E).values())
          if (m.getRank()==1) tops.add(getCombinedMatch(m, E, dirt));
//    }
    return tops;
  }
  // locate the poorest performing matching for each engine to server as the LOD
  protected Map<Engine, PeptideMatch> getEngineDirt()
  {
    Map<Engine, PeptideMatch> engine_dirt = new HashMap<>();
    for (Engine E : mPeptideEngineMatch.columnKeySet())
    {
      PeptideMatch dirt=null; double multiple=0.5;
      if (mPeptideEngineMatch.column(E).values().size()>1)
      {
        for (PeptideMatch d : mPeptideEngineMatch.column(E).values())
          if (dirt==null || dirt.getRank()<d.getRank()) dirt=d;
      }
      else
      {
        // the engine outputs only the top match for the scan. Have to fake one!
        dirt = Tools.front(mPeptideEngineMatch.column(E).values()); multiple=0.01;
      }

      if (dirt!=null)
      {
        PeptideMatch d = PSMs.clone2LOD(dirt, multiple);
        d.setRejected(true);
        engine_dirt.put(E, d);
      }
    }

    return engine_dirt;
  }
  protected PeptideMatch getCombinedMatch(String peptide, Map<Engine, PeptideMatch> engine_dirt)
  {
    PeptideMatch top = PSMs.clone(Tools.front(mPeptideEngineMatch.row(peptide).values()));

    // merge the scores
    top.getScoreMap().clear(); double voted=0d;
    for (Engine E : mPeptideEngineMatch.columnKeySet())
    {
      PeptideMatch m = mPeptideEngineMatch.get(peptide, E);

      if (m!=null && !m.isRejected()) voted+=1d;
      else
      {
        m = engine_dirt.get(E); m.setRejected(true);
      }

      // signal that it has already appear in the output
      m.setRank(m.getRank()*-1);
      for (String scr : m.getScoreMap().keySet())
        if (scr.charAt(0)!='^') top.addScore(scr, m.isRejected()?(Stats.d2d(m.getScore(scr),2)):m.getScore(scr));
    }
    top.addScore("Vote", voted); top.setRank(1);

    return top;
  }
  // fill in the missing scores
  protected PeptideMatch getCombinedMatch(PeptideMatch match, Engine E, Map<Engine, PeptideMatch> engine_dirt)
  {
    PeptideMatch top = PSMs.clone(match);

    // set the vote = 1
    top.addScore("Vote", 1); top.setRank(1);

    // merge the scores
    for (Engine engine : engine_dirt.keySet())
      if (!engine.equals(E))
        for (String scr : engine_dirt.get(engine).getScoreMap().keySet())
          // truncate the score to 1 deci point to indicate this is a fake one
          if (scr.charAt(0)!='^') top.addScore(scr, Stats.d2d(engine_dirt.get(engine).getScore(scr), 2));

    return top;
  }

  // Getters and Setters
  public SpectrumIdentifier getEffectiveSpectrumID()
  {
    // TODO need to create a composite spec ID if we are dealing with super consensus of multiple related spectra such as charge variants
    return mRunScanSpectrum!=null?Tools.front(mRunScanSpectrum.values()):null;
  }
  public Table<String, String, SpectrumIdentifier> getSpectrumIDs()        { return mRunScanSpectrum; }
  public Table<String, Engine, PeptideMatch>       getPeptideEngineMatch() { return mPeptideEngineMatch; }

  /** Go thro the alignments and rank them by some kind of criteria
   *
   * @return the object itself
   */
  public MetaPeptideMatch rank()
  {
    if (Tools.isSet(mPeptideEngineMatch))
      for (String peptide : getPeptideEngineMatch().rowKeySet())
      {
        double r = 0;
        for (PeptideMatch m : getPeptideEngineMatch().row(peptide).values()) r+=m.getRank();
        r/=getPeptideEngineMatch().row(peptide).size();
        mVoteRankSeqIndex.put(getPeptideEngineMatch().row(peptide).size()*-1, r, peptide);
      }

    return this;
  }
  // add the PSM one at a time
  public MetaPeptideMatch addAll(SpectrumIdentifier ms, Collection<PeptideMatch> mm, Engine search)
  {
    if (Tools.isSet(mm))
      for (PeptideMatch m : mm) add(ms, m, search);

    return this;
  }
  public MetaPeptideMatch add(SpectrumIdentifier ms, PeptideMatch m, Engine engine)
  {
    if (ms!=null && m!=null && engine!=null)
    {
      mRunScanSpectrum.put(Strs.stripLastOf(ms.getSpectrumFile().get(), '.'), LCMSMS.toScanStr(ms.getScanNumbers()), ms);
//      String       p = PSMs.toNumModSequence(m);
      String       p = m.toBarePeptide().toSymbolString();

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
  public static Map<String, MetaPeptideMatch> add(Map<String, MetaPeptideMatch> row_matches, Engine engine,
                                            Multimap<SpectrumIdentifier, PeptideMatch> id_match)
  {
    if (row_matches==null) row_matches = new TreeMap<>();

    if (engine!=null && Tools.isSet(id_match))
      for (SpectrumIdentifier id : id_match.keySet())
      {
        String   row = id.getSpectrumFile().get()+"#"+ LCMSMS.toScanStr(id.getScanNumbers());
        MetaPeptideMatch E = row_matches.get(row);
        if (E==null) { E = new MetaPeptideMatch(); row_matches.put(row, E); }
        // add the contents of the matches
        for (PeptideMatch m : id_match.get(id))
        {
          E.add(id, m, engine);
        }
      }

    return row_matches;
  }
  @Override
  public String toString()
  {
    StringBuffer buf = new StringBuffer();

    buf.append(getEffectiveSpectrumID().getSpectrum()+"\n");
    if (Tools.isSet(getPeptideEngineMatch()))
    {
      buf.append("Peptide");
      for (Engine engine : getPeptideEngineMatch().columnKeySet()) buf.append("\t"+engine.getName());
      buf.append("\n");
      for (String peptide : getPeptideEngineMatch().rowKeySet())
      {
        for (Engine engine : getPeptideEngineMatch().columnKeySet())
        {
          PeptideMatch m = getPeptideEngineMatch().get(peptide, engine);
          buf.append((m!=null?m.getRank()+":"+Tools.d2s(m.getScore(engine.getCanonicalScore()), 2):"-")+"\t");
        }
        buf.append(peptide+"\n");
      }
    }
    return buf.toString();
  }
  public static void rank(Collection<MetaPeptideMatch> matches)
  {
    for (MetaPeptideMatch m : matches) m.rank();
  }
  // generate a pooled series of worst performing matches
  public static RowSortedTable<Double, Engine, PeptideMatch> getEngineDirts(Map<SpectrumIdentifier, MetaPeptideMatch> meta, double pool)
  {
    RowSortedTable<Double, Engine, PeptideMatch> rt_engine_dirt = TreeBasedTable.create();
    for (SpectrumIdentifier id : meta.keySet())
    {
      Map<Engine, PeptideMatch> dirt = meta.get(id).getEngineDirt();
      if (Tools.isSet(dirt))
        for (Engine E : dirt.keySet())
          if (dirt.containsKey(E))
          {
            Double pos = (double )Math.round(id.getScanNumbers().getFirst().getValue() / pool);
            if (!rt_engine_dirt.contains(pos, E) ||
               E.compareTo(rt_engine_dirt.get(pos, E).getScoreMap(), dirt.get(E).getScoreMap())<0)
              rt_engine_dirt.put(pos, E, dirt.get(E));
          }
    }
    return rt_engine_dirt;
  }
  public static Multimap<SpectrumIdentifier, PeptideMatch> getTopPeptideMatches(
    Map<SpectrumIdentifier, MetaPeptideMatch> meta, boolean keepConsensus)
  {
    double pool = 100d;
    RowSortedTable<Double, Engine, PeptideMatch> rt_engine_dirt = getEngineDirts(meta, pool);
//    for (SpectrumIdentifier id : meta.keySet())
//    {
//      Map<Engine, PeptideMatch> dirt = meta.get(id).getEngineDirt();
//      if (Tools.isSet(dirt))
//        for (Engine E : dirt.keySet())
//          rt_engine_dirt.put(id.getScanNumbers().getFirst().getValue(), E, dirt.get(E));
//    }

    Multimap<SpectrumIdentifier, PeptideMatch> tops = HashMultimap.create();
    // loop thro the matches
    Double last_pos=null, pos=null;
    Map<Engine, PeptideMatch> dirt=null;
    for (SpectrumIdentifier id : meta.keySet())
    {
      pos = (double )Math.round(id.getScanNumbers().getFirst().getValue()/pool);
      if (last_pos==null || pos!=last_pos)
      {
        dirt     = interpolate(rt_engine_dirt, pos);
        last_pos = pos;
      }

//      if (id.getScanNumbers().getFirst().getValue()==1914)
//      {
//        System.out.println();
//      }
      if (keepConsensus) tops.putAll(id, meta.get(id).getTopPeptideMatches(dirt));
      else               tops.put(   id, meta.get(id).getTopPeptideMatch(  dirt));
    }

    return tops;
  }
  public static Map<Engine, PeptideMatch> interpolate(RowSortedTable<Double, Engine, PeptideMatch> dirts, double rt)
  {
    Map<Engine, PeptideMatch> dirt = new HashMap<>();

    Table<Double, Engine, PeptideMatch> sub = Tools.interpolate(dirts, rt);
    if (sub!=null)
    {
      // looking for the worst candidate
      for (Engine E : sub.columnKeySet())
        if (Tools.isSet(sub.column(E)))
          dirt.put(E, PSMs.worst(E, sub.column(E).values()));
    }

    return dirt;
  }
}
