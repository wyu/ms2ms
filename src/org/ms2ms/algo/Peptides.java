package org.ms2ms.algo;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.expasy.mzjava.core.mol.SymbolSequence;
import org.expasy.mzjava.proteomics.mol.Protein;
import org.expasy.mzjava.proteomics.mol.modification.ModAttachment;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.util.Collection;
import java.util.List;

/**
 * Created by yuw on 9/14/2015.
 */
public class Peptides
{
  public static String toNumModSequence(PeptideMatch p)
  {
    if (p==null) return null;

    final StringBuilder sb = new StringBuilder();

    sb.append((Tools.isSet(p.getModifications(ModAttachment.nTermSet)) ||
               Tools.isSet(p.getModifications(0, ModAttachment.sideChainSet))) ?
               p.getSymbol(0).getSymbol().toLowerCase() :
               p.getSymbol(0).getSymbol().toUpperCase());

    for (int i = 1; i < p.size()-1; i++)
    {
      sb.append(Tools.isSet(p.getModifications(i, ModAttachment.sideChainSet)) ?
                p.getSymbol(i).getSymbol().toLowerCase() :
                p.getSymbol(i).getSymbol().toUpperCase());
    }

    sb.append((Tools.isSet(p.getModifications(ModAttachment.cTermSet)) ||
               Tools.isSet(p.getModifications(p.size()-1, ModAttachment.sideChainSet))) ?
               p.getSymbol(p.size()-1).getSymbol().toLowerCase() :
               p.getSymbol(p.size()-1).getSymbol().toUpperCase());

    return sb.toString();
  }
  // +229.163HMK+229.163K+229.163HAK+229.163K+229.163MK+229.163K+229.163QMK+229.163K+229.163
  public static PeptideMatch fromNumModSequence(String p)
  {
    if (!Strs.isSet(p)) return null;

    List<String> items = Strs.splits(p, "[+-.\\d]+");
    // an array of: "", +229.163, HMK, +229.163, K, +229.163, HAK, +229.163, K, +229.163, MK, +229.163, K, +229.163, QMK, +229.163, K, +229.163
    String backbone = null;
    if (Tools.isSet(items))
    {
      // merge to the backbone first
      for (int i=0; i<items.size(); i+=2)
        backbone = Strs.extend(backbone, items.get(i), "");

      PeptideMatch m = new PeptideMatch(backbone); backbone=null;
      for (int i=0; i<items.size(); i+=2)
      {
        backbone = Strs.extend(backbone, items.get(i), "");
        // update the peptide match

        // the mod string
        if (!Strs.isSet(backbone))
        {
          // N-term mod
          m.addModificationMatch(ModAttachment.N_TERM, Double.valueOf(items.get(i+1)));
        }
        else if (i+1<items.size())
        {
          // Side-chain mod
          m.addModificationMatch(backbone.length()-1, Double.valueOf(items.get(i+1)));
        }
      }
      return m;
    }
    return null;
  }
  public static PeptideMatch addScore(PeptideMatch m, String t, Double s)
  {
    if (m!=null && s!=null && Strs.isSet(t)) m.addScore(t, s);
    return m;
  }
  public static PeptideMatch addScore(PeptideMatch m, String t, Integer s)
  {
    if (m!=null && s!=null && Strs.isSet(t)) m.addScore(t, s);
    return m;
  }
  public static <P extends SymbolSequence> Multimap<String, P> toSequenceMap(Collection<P>... proteins)
  {
    if (!Tools.isSet(proteins)) return null;
    // create a map
    Multimap<String, P> seq_prots = HashMultimap.create();
    for (Collection<P> ps : proteins)
      for (P p : ps) seq_prots.put(p.toSymbolString(), p);

    return seq_prots;
  }
  public static boolean hasScoreAt(PeptideMatch m, String tag, double score)
  {
    // anywhere among the ranks
    return (m!=null && m.hasScore(tag) && m.getScore(tag)==score);
  }
  public static boolean hasScoreAtBelow(PeptideMatch m, String tag, double score)
  {
    // anywhere among the ranks
    return (m!=null && m.hasScore(tag) && m.getScore(tag)<=score);
  }
}
