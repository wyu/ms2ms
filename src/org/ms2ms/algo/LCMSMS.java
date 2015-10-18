package org.ms2ms.algo;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.TreeBasedTable;
import com.google.common.collect.TreeMultimap;
import org.expasy.mzjava.core.ms.spectrum.ScanNumber;
import org.expasy.mzjava.core.ms.spectrum.ScanNumberList;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatchComparator;
import org.expasy.mzjava.proteomics.ms.ident.PeptideProteinMatch;
import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import org.ms2ms.math.Stats;
import org.ms2ms.r.Dataframe;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentification;

import java.util.*;

/**
 * Created by yuw on 9/20/2015.
 */
public class LCMSMS
{
  public static String toScanStr(ScanNumberList scans)
  {
    String scan=null;
    if (Tools.isSet(scans))
    {
      for (ScanNumber s : scans)
        scan = Strs.extend(scan, s.getValue()+"", "+");
    }
//    if (Strs.equals(scan, "-1"))
//    {
//      System.out.println();
//    }
    return scan;
  }
  public static String toRunScan(SpectrumIdentifier id)
  {
    if (id==null) return null;

    String scan=null;
    if (id.getSpectrumFile().isPresent()) Strs.extend(scan, id.getSpectrumFile().get(),"");

    return Strs.extend(scan, toScanStr(id.getScanNumbers()), "#");
  }
  public static Multimap<SpectrumIdentifier, PeptideMatch> byRun(Multimap<SpectrumIdentifier, PeptideMatch> id_match, String... runs)
  {
    if (Tools.isSet(id_match) && Tools.isSet(runs))
    {
      Iterator<SpectrumIdentifier> itr = id_match.keySet().iterator();
      while (itr.hasNext())
      {
        if (!Strs.isA(Strs.stripLastOf(itr.next().getSpectrumFile().get(), '.'), runs)) itr.remove();
      }
    }
    return id_match;
  }
  public static Multimap<SpectrumIdentifier, PeptideMatch> byQ(Multimap<SpectrumIdentifier, PeptideMatch> id_match, String score, double q)
  {
    if (Tools.isSet(id_match))
    {
      Iterator<SpectrumIdentifier> itr = id_match.keySet().iterator();
      while (itr.hasNext())
      {
        boolean qualified = false;
        for (PeptideMatch m : id_match.get(itr.next()))
        {
          if (m.getScore(score)<=q) { qualified=true; break; }
        }
        if (!qualified) itr.remove();
      }
    }
    return id_match;
  }
  // consider only the top-ranked
  public static <K> Dataframe cut(Multimap<K, PeptideMatch> id_match, String score, double... qvals)
  {
    if (Tools.isSet(id_match))
    {
      id_match = byRank(id_match, 1);
      Dataframe d = new Dataframe("PSMs thresholded by " + score);
      int row=0;
      for (double q : qvals)
      {
        System.out.print("Considering the threshold. ");
        long hits=0, decoys=0, inProtein=0, inProteinTop=0, inProteinTopQual=0;
        for (K id : id_match.keySet())
        {
          boolean inP=false;
          for (PeptideMatch m : id_match.get(id))
          {
            // anywhere among the ranks
            boolean inProt=Peptides.hasScoreAt(m, "inProtein", 1);
            if (inProt) inP=true;
            // only count the top ranked
            if (m.getRank()==1 && inProt)
            {
              inProteinTop++;
              if (m.getScore(score)<=q) inProteinTopQual++;
            }
            if (m.getRank()==1 && m.hasScore(score) && m.getScore(score) <= q)
            {
              if (m.getProteinMatches().get(0).getHitType().equals(PeptideProteinMatch.HitType.DECOY)) decoys++;
              else hits++;
//              break;
            }
          }
          if (inP) inProtein++;
        }
        // populate the data frame
        d.put(row, "Threshold", q);
        d.put(row, "Hits", hits);
        d.put(row, "Decoys", decoys);
        d.put(row, "inProtein", inProtein);
        d.put(row, "inProteinTop", inProteinTop);
        d.put(row, "inProteinTopQual", inProteinTopQual);
        d.put(row, "MSMS", id_match.keySet().size());

        System.out.println("Threshold: " + q + ", Hits: " + hits + ", Decoys: " + decoys + ", MSMS: " + id_match.keySet().size());
        row++;
      }
      return d;
    }
    return null;
  }
  public static Multimap<SpectrumIdentifier, PeptideMatch> rank(Multimap<SpectrumIdentifier, PeptideMatch> id_match, String score)
  {
    if (Tools.isSet(id_match))
    {
      // sort out the rank of the matches first
      TreeMultimap<Double, PeptideMatch> score_match = TreeMultimap.create();
      for (SpectrumIdentifier id : id_match.keySet())
      {
        List<PeptideMatch> matches = new ArrayList<PeptideMatch>(id_match.get(id));
        Collections.sort(matches, new PeptideMatchComparator(score));
        for (int i=0; i< matches.size(); i++) matches.get(i).setRank(i+1);
/*
        score_match.clear();
        for (PeptideMatch m : id_match.get(id)) score_match.put(m.getScore(score)*-1d, m);
        int rank=1;
        for (Double scr : score_match.keySet())
        {
          for (PeptideMatch m : score_match.get(scr)) m.setRank(rank);
          rank++;
        }
*/
      }
    }
    return id_match;
  }
  public static <K> Multimap<K, PeptideMatch> byRank(Multimap<K, PeptideMatch> id_match, int n)
  {
    if (Tools.isSet(id_match))
    {
      Multimap<K, PeptideMatch> remain = HashMultimap.create();
      for (K id : id_match.keySet())
      {
        for (PeptideMatch m : id_match.get(id))
          if (m.getRank()<=n) remain.put(id, m);
      }
      return remain;
    }
    return id_match;
  }
  public static Multimap<String, PeptideMatch> toScanMatch(Multimap<SpectrumIdentifier, PeptideMatch> id_match)
  {
    if (id_match==null) return null;

    Multimap<String, PeptideMatch> scan_match = TreeMultimap.create();
    for (SpectrumIdentifier id : id_match.keySet())
      scan_match.putAll(LCMSMS.toScanStr(id.getScanNumbers()), id_match.get(id));

    return scan_match;
  }
  public static Double parseNominalPrecursorMz(String s)
  {
    if (Strs.isSet(s))
    {
      String[] items = Strs.split(s.substring(s.indexOf("ms"), s.indexOf('[')), ' ', true);
      // grab the MS2 info and change the precursor to that of the MS2
      return Stats.toDouble(Strs.split(items[1], '@', true)[0]);
    }
    return null;
  }
  // check the protein sequences to see if any peptide is a match
  public static Multimap<String, PeptideMatch> byProteinSequence(Collection<String> sequences, Collection<PeptideMatch>... hits)
  {
    if (!Tools.isSet(hits) || !Tools.isSet(hits[0])) return null;

    // let's build the index first
    long inP=0, totals=0;
    Multimap<String, PeptideMatch> seq_match = Peptides.toSequenceMap(hits);
    for (String sequence : sequences)
      for (String peptide : seq_match.keySet())
      {
        totals+=seq_match.get(peptide).size();
        if (sequence.indexOf(peptide)>= 0)
          for (PeptideMatch m : seq_match.get(peptide)) { inP++; m.addScore("inProtein", 1); }
      }

    System.out.println("Matches (inProtein/Totals): " + inP+"/"+totals);
    return seq_match;
  }
  // merge the PSMs from Bs onto As, by their precursors information. The scan # was missing from the decoy PSMs export from SequestHT/PD
/*
  public static Multimap<String, PeptideMatch> mergeByPrecursor(Multimap<String, PeptideMatch> As, Multimap<String, PeptideMatch> Bs)
  {
    if (!Tools.isSet(As)) return Bs;
    if (!Tools.isSet(Bs)) return As;

    // let's create a index first
    Table<Double, Integer, PeptideMatch> mass_match = TreeBasedTable.create();
    for (PeptideMatch m : Bs.values())
      mass_match.put(m.getNeutralPeptideMass(), m);

    for (String id : As.keySet())
    {
      PeptideMatch m1 = Tools.front(As.get(id));
      Collection<PeptideMatch> found = mass_match.get(m1.getNeutralPeptideMass());
    }
    return null;
  }
*/
}
