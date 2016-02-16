package org.ms2ms.algo;

import com.google.common.base.*;
import com.google.common.collect.*;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import org.apache.commons.collections.map.HashedMap;
import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.expasy.mzjava.proteomics.mol.AminoAcid;
import org.expasy.mzjava.proteomics.mol.Peptide;
import org.expasy.mzjava.proteomics.mol.modification.ModAttachment;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.mol.modification.unimod.UnimodManager;
import org.expasy.mzjava.proteomics.ms.fragment.PeptideFragmentAnnotator;
import org.expasy.mzjava.proteomics.ms.fragment.PeptideFragmenter;
import org.expasy.mzjava.proteomics.ms.ident.*;
import org.expasy.mzjava.proteomics.ms.spectrum.PepFragAnnotation;
import org.expasy.mzjava.proteomics.ms.spectrum.PeptideSpectrum;
import org.ms2ms.data.ms.Engine;
import org.ms2ms.data.ms.MetaPeptideMatch;
import org.ms2ms.io.PsmReaders;
import org.ms2ms.mzjava.NumModMatchResolver;
import org.ms2ms.mzjava.SpectrumIdentifierByRunScan;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.TabFile;
import org.ms2ms.utils.Tools;

import javax.annotation.Nonnull;
import java.io.DataOutput;
import java.io.IOException;
import java.util.*;

/** Routines about the peptide-spectrum-matches.
 *
 * Created by yuw on 10/19/2015.
 */
public class PSMs
{
  public static final String SCR_CANONICAL = "Canonical Score";
  public static final String SCR_DELTA     = "Delta";

  // the comparators
  public static class DesendCanonicalScorePeptideMatch implements Comparator<PeptideMatch>
  {
    String score = SCR_CANONICAL;
    public DesendCanonicalScorePeptideMatch()         { super(); };
    public DesendCanonicalScorePeptideMatch(String s) { score = s; }
    public int compare(PeptideMatch o1, PeptideMatch o2)
    {
      return o1!=null && o2!=null ? Double.compare(o2.getScore(score), o1.getScore(score)):0;
    }
  }
  public static NumModMatchResolver sNumModResolver = new NumModMatchResolver();
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
        if (Strs.isSet(items.get(i)))
          for (char c : items.get(i).toCharArray())
            if (Peptides.sAA.indexOf(c)>=0)
              backbone = Strs.extend(backbone, c+"", "");

      PeptideMatch m = new PeptideMatch(backbone); backbone=null;
      for (int i=0; i<items.size(); i+=2)
      {
        if (Strs.isSet(items.get(i)))
          for (char c : items.get(i).toCharArray())
            if (Peptides.sAA.indexOf(c)>=0)
              backbone = Strs.extend(backbone, c+"", "");
//        backbone = Strs.extend(backbone, items.get(i), "");
        // update the peptide match

        // the mod string
        if (!Strs.isSet(backbone))
        {
          try
          {
            List<String> nmods = Strs.split(items.get(i+1), "(?=[(+)(\\-)])", true);
            // N-term mod
            for (String nmod : nmods)
              m.addModificationMatch(ModAttachment.N_TERM, Double.valueOf(nmod));
          }
          catch (Exception e)
          {
            e.printStackTrace();
          }
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

  /** FragmentAnnotator is a facade that highly simplifies and hides the complexity of annotating a PeakList given a peptide.
   *  Annotating a peak list requires a theoretical spectra to be generated, the peaks of the theoretical and query spectra
   *  to be aligned and reporting annotations back to the PeakList. Here is a java snippet to annotate a PeakList from a Peptide:
   *  http://mzjava.expasy.org/daily/wiki/Recipe-4.1.html
   *
   * @param spectrum
   * @param peptide
   * @return
   */
  public static PeakList<PepFragAnnotation> annotate(PeakList<PepFragAnnotation> spectrum, Peptide peptide)
  {
    // the fragmenter is responsible for generating a theoretical spectrum from the peptide
    PeptideFragmenter fragmenter =
        new PeptideFragmenter(EnumSet.of(IonType.b, IonType.y), PeakList.Precision.DOUBLE);

    // the annotator needs to delegate to the fragmenter the fragmentation process and to the internal aligner the tolerance for aligning peaks
    PeptideFragmentAnnotator annotator = new PeptideFragmentAnnotator(fragmenter, new AbsoluteTolerance(0.1));

    // the peptide as a source for theoretical fragments with annotations
    // Peptide peptide = Peptide.parse("QVHPDTGISSK");

    // start process: fragmenting peptide, aligning to the spectrum and reporting matched ions back to the spectrum
    return annotator.annotate(spectrum, peptide, 3);
  }
  public static PeptideSpectrum annotate(PeptideFragmentAnnotator annotator, MsnSpectrum spectrum, PeptideMatch match)
  {
    // convert the spectrum to the suitable format for the annotator
    Peptide      peptide = match.toPeptide(sNumModResolver);
    PeptideSpectrum msms = new PeptideSpectrum(peptide, spectrum.getPrecursor().getCharge());
    for (int i=0; i<spectrum.size(); i++)
      msms.add(spectrum.getMz(i), spectrum.getIntensity(i));

    // start process: fragmenting peptide, aligning to the spectrum and reporting matched ions back to the spectrum
    return annotator.annotate(msms, peptide, spectrum.getPrecursor().getCharge());
  }
  public static boolean hasFragType(Collection<PepFragAnnotation> as, IonType... types)
  {
    if (Tools.isSet(as))
      for (PepFragAnnotation A : as)
        for (IonType type : types)
          if (type.equals(A.getIonType())) return true;

    return false;
  }

  /** Merge the PSMs by the scan.
   *
   * @param psms are the search outputs from two or more engines
   */
  public static Map<SpectrumIdentifier, MetaPeptideMatch> merge(
    Map<SpectrumIdentifier, MetaPeptideMatch>  meta, Engine search,
    Multimap<SpectrumIdentifier, PeptideMatch> psms)
  {
    if (meta==null) meta = new HashMap<>();
    // create an index
    Map<String, SpectrumIdentifier> scan_id = new HashMap<>();
    for (SpectrumIdentifier id : meta.keySet()) scan_id.put(id.getSpectrum(), id);

    // loop thro the PSMs
    for (SpectrumIdentifier id : psms.keySet())
    {
      MetaPeptideMatch mm = meta.get(scan_id.get(id.getSpectrum()));
      if (mm==null)
      {
        mm = new MetaPeptideMatch();
        scan_id.put(id.getSpectrum(), id); meta.put(id, mm);
      }
      // fill in the matches
      mm.addAll(id, psms.get(id), search);
    }

    return meta;
  }
  public static Multimap<SpectrumIdentifier, PeptideMatch> accumulate(
      Multimap<SpectrumIdentifier, PeptideMatch> As,
      Multimap<SpectrumIdentifier, PeptideMatch> Bs, int tops, String variety, String score)
  {
    if (!Tools.isSet(Bs)) return As;
    if (As==null) As = HashMultimap.create();

    System.out.print(Bs.size() + " total PSMs from " + Bs.keySet().size() + "MS/MS");
    // combine the PSMs
    PSMs.DesendCanonicalScorePeptideMatch sorter = new PSMs.DesendCanonicalScorePeptideMatch(score);
    List<PeptideMatch>                   matches = new ArrayList<>();
    for (SpectrumIdentifier id : Bs.keySet())
    {
      // System.out.println(" --> " + As.size() + " total PSMs from " + As.keySet().size() + "MS/MS");
      // tag the PSM to the search condition
      for (PeptideMatch m : Bs.get(id)) m.addScore(variety, 1);
      // pool the matches for the scan together
      matches.clear(); matches.addAll(Bs.get(id));
      if (As.containsKey(id)) matches.addAll(As.get(id));
      // rank the matches by the canonical score in desending order
      Collections.sort(matches, sorter);
      // keep only the top few to preserve the memory footprint
      As.removeAll(id);
      As.putAll(id, matches.subList(0, tops<=matches.size()?tops:matches.size()));
    }
    System.out.println(" --> " + As.size() + " total PSMs from " + As.keySet().size() + " MS/MS");

    return As;
  }
  public static PeptideMatch clone(@Nonnull PeptideMatch m)
  {
    PeptideMatch M = new PeptideMatch(m.toSymbolString());

    // copy the contents
    M.setRank(m.getRank());
    M.setNumMatchedIons(m.getNumMatchedIons());
    M.setNumMissedCleavages(m.getNumMissedCleavages());
    M.setTotalNumIons(m.getTotalNumIons());
    M.setNeutralPeptideMass(m.getNeutralPeptideMass());
    M.setMassDiff(m.getMassDiff());
    M.setRejected(m.isRejected());
    M.addProteinMatches(m.getProteinMatches());

    Collection<ModificationMatch> all = m.getModifications(ModAttachment.all);
    if (Tools.isSet(all))
      for (ModificationMatch mod : all)
        if (ModAttachment.SIDE_CHAIN.equals(mod.getModAttachment())) M.addModificationMatch(mod.getPosition(), mod);
        else                                                         M.addModificationMatch(mod.getModAttachment(), mod);

    if (m.getScoreMap()!=null)
      for (String score : m.getScoreMap().keySet())
        M.addScore(score, m.getScore(score));

    return M;
  }

  /** Create a clone of the match and set its scores to the multiple of the original
   *
   * @param m is the original match to be cloned
   * @param multiple is the score multiple required
   * @return the cloned match with the given score multiple
   */
  public static PeptideMatch clone2LOD(@Nonnull PeptideMatch m, double multiple)
  {
    // look for the lowest rank for this engine
    PeptideMatch dummy = PSMs.clone(m); dummy.getScoreMap().clear();
    // set the score multiples
    for (String score : m.getScoreMap().keySet())
      dummy.addScore(score, m.getScore(score)*multiple);

    return dummy;
  }
  public static PeptideMatch worst(@Nonnull Engine E, @Nonnull Collection<PeptideMatch> matches)
  {
    PeptideMatch worst=null;
    for (PeptideMatch match : matches)
      if (worst==null || (E.compareTo(worst.getScoreMap(), match.getScoreMap())<0)) worst=match;

    return worst;
  }

  public static Multimap<String, PeptideMatch> toBackboneMatchMap(@Nonnull Collection<PeptideMatch> matches)
  {
    if (matches==null) return null;

    Multimap<String, PeptideMatch> seq_match = HashMultimap.create();
    for (PeptideMatch match : matches)
      seq_match.put(match.toSymbolString(), match);

    return seq_match;
  }
/*
  */
/**@param mm are the input peptide matches
   * @param multiple is the score multiple required
   * @return a cloned match where the worst of the
   *//*

  public static PeptideMatch worst2LOD(@Nonnull Collection<PeptideMatch> mm, double multiple)
  {
    // look for the lowest rank for this engine
    PeptideMatch dummy = PSMs.clone(Tools.front(mm)); dummy.getScoreMap().clear();
    // looking for the worst
    Map<String, Double> bottom_feeder = new HashMap<>();
    for (PeptideMatch m : mm)
      for (String score : m.getScoreMap().keySet())
        if (!bottom_feeder.containsKey(score) || bottom_feeder.get(score)<m.getScore(score))
          bottom_feeder.put(score, m.getScore(score));

    // set the score multiples
    for (String score : bottom_feeder.keySet())
      dummy.addScore(score, bottom_feeder.get(score)*multiple);

    return dummy;
  }
*/
}
