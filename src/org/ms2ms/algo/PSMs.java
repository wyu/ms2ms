package org.ms2ms.algo;

import com.google.common.collect.*;
import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.expasy.mzjava.proteomics.mol.Peptide;
import org.expasy.mzjava.proteomics.mol.modification.ModAttachment;
import org.expasy.mzjava.proteomics.ms.fragment.PeptideFragmentAnnotator;
import org.expasy.mzjava.proteomics.ms.fragment.PeptideFragmenter;
import org.expasy.mzjava.proteomics.ms.ident.*;
import org.expasy.mzjava.proteomics.ms.spectrum.PepFragAnnotation;
import org.expasy.mzjava.proteomics.ms.spectrum.PeptideSpectrum;
import org.ms2ms.data.ms.Engine;
import org.ms2ms.data.ms.MetaPeptideMatch;
import org.ms2ms.data.ms.PSM;
import org.ms2ms.math.Stats;
import org.ms2ms.mzjava.NumModMatchResolver;
import org.ms2ms.r.Dataframe;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.util.*;

/** Routines about the key-spectrum-matches.
 *
 * Created by yuw on 10/19/2015.
 */
public class PSMs
{
  public static final String SCR_CANONICAL = "CanonicalScore";
  public static final String SCR_DELTA     = "DeltaScore";

  public static Map<String, PeptideMatch> mods_match = new HashMap<>();

  // the comparators
  public static class DesendScorePeptideMatch implements Comparator<PeptideMatch>
  {
    String score = SCR_CANONICAL;
    public DesendScorePeptideMatch()         { super(); };
    public DesendScorePeptideMatch(String s) { score = s; }
    public int compare(PeptideMatch o1, PeptideMatch o2)
    {
      if (!o1.hasScore(score) || !o2.hasScore(score))
      {
//        System.out.println();
      }
      return o1!=null && o2!=null ? Double.compare(o2.getScore(score), o1.getScore(score)):0;
    }
  }

  public static void verify(PeptideMatch m, MsnSpectrum ms2, MsnSpectrum ms3, PeptideFragmenter fragmenter)
  {
//    PeptideFragmenter fragmenter = new PeptideFragmenter(EnumSet.of(IonType.b, IonType.y), PeakList.Precision.DOUBLE);
    Peptide peptide = Peptide.parse("GYDSPR");

// generating key b and y fragments with charges +1 and +2
    PeptideSpectrum peptideSpectrum = fragmenter.fragment(m.toPeptide(), 2);

  }
  public static Dataframe alignment(Dataframe consensus, Engine engine, Multimap<SpectrumIdentifier, PeptideMatch> matches)
  {
    return alignment(consensus, engine, engine.getName(), matches);
  }
  // going thro the outputs from each engine and produce the alignment by run-scan-backbone
  public static Dataframe alignment(Dataframe consensus, Engine engine, String name, Multimap<SpectrumIdentifier, PeptideMatch> matches)
  {
    if (!Tools.isSet(matches)) return consensus;
    if (consensus==null) consensus = new Dataframe();

    // go over each of the key matches
    for (SpectrumIdentifier id : matches.keySet())
    {
      if (!id.getName().isPresent())
      {
        throw new RuntimeException("Scan name not specified!");
      }
      String run = Strs.split(id.getName().get(), '#')[0];
      for (PeptideMatch m : matches.get(id))
      {
        String row = id.getName().get()+"^"+m.toSymbolString();
        // deposit some of the property cols
        consensus.put(row, "Run",     run);
        consensus.put(row, "Scan",    id.getScanNumbers().getFirst().getValue());
        consensus.put(row, "Sequence", m.toSymbolString());
        consensus.put(row, "Decoy",    Tools.equals(PeptideProteinMatch.HitType.DECOY, Tools.front(m.getProteinMatches()).getHitType()));

        if (id.getPrecursorMz(  ).isPresent()) consensus.put(row, "m/z",id.getPrecursorMz().get());
        if (id.getAssumedCharge().isPresent()) consensus.put(row, "z",  id.getAssumedCharge().get());

        Integer n = consensus.cell(row, "Mods")!=null && consensus.cell(row, "Mods") instanceof Integer ? (Integer )consensus.cell(row, "Mods") : null;
        if (m.getModificationCount()>0 && (n==null || n<m.getModificationCount())) {
          consensus.put("ModN", m.getModificationCount());
          consensus.put("Peptide", PSMs.toNumModSequence(m));
        }

        consensus.put(row, name,       m.getScore(engine.getCanonicalScore()));
      }
    }
    System.out.println("--> " + consensus.size());

    return consensus;
  }
  // going thro the outputs from each engine and produce the alignment by run-scan-backbone
  public static Table<SpectrumIdentifier, Engine, PeptideMatch> alignment(
      Table<SpectrumIdentifier, Engine, PeptideMatch> consensus, Engine engine, Multimap<SpectrumIdentifier, PeptideMatch> matches)
  {
    if (!Tools.isSet(matches)) return consensus;
    if (consensus==null) consensus = HashBasedTable.create();

    Map<String, SpectrumIdentifier> runscan2id = new HashMap<>();
    if (Tools.isSet(consensus))
      for (SpectrumIdentifier id : consensus.rowKeySet())
        runscan2id.put(id.getSpectrum(), id);

      // go over each of the key matches
    for (SpectrumIdentifier id : matches.keySet())
    {
      if (!id.getName().isPresent())
      {
        throw new RuntimeException("Scan name not specified!");
      }
      for (PeptideMatch m : matches.get(id))
      {
        // deposit some of the property cols
        SpectrumIdentifier id2 = runscan2id.get(id.getSpectrum());
        if (id2==null)
        {
          runscan2id.put(id.getSpectrum(), id);
        }
        consensus.put(runscan2id.get(id.getSpectrum()), engine, m);
      }
    }
    System.out.println("--> " + consensus.size());

    return consensus;
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
  public static PeptideMatch fromNumModSequence(String p)
  {
    if (!Strs.isSet(p)) return null;

//    // always have much fewer distinct peptides than PSMs
//    if (mods_match.containsKey(p)) return mods_match.get(p);

    List<String> items = Strs.splits(p, Strs.PTN_SIGNS_DGT);
    // an array of: "", +229.163, HMK, +229.163, K, +229.163, HAK, +229.163, K, +229.163, MK, +229.163, K, +229.163, QMK, +229.163, K, +229.163
    // use StringBuffer instead of String to avoid GC overhead problem? WYU, 20160216
    StringBuffer backbone = new StringBuffer(48);
    if (Tools.isSet(items))
    {
      // merge to the backbone first
      for (int i=0; i<items.size(); i+=2)
        if (Strs.isSet(items.get(i)))
          for (char c : items.get(i).toCharArray())
            if (Peptides.sAA.indexOf(c)>=0)
              backbone.append(c);

      PeptideMatch m = new PeptideMatch(backbone.toString()); backbone=new StringBuffer(48);
      for (int i=0; i<items.size(); i+=2)
      {
        if (Strs.isSet(items.get(i)))
          for (char c : items.get(i).toCharArray())
            if (Peptides.sAA.indexOf(c)>=0)
              backbone.append(c);

        // the mod string
        try
        {
          // hope the pre-compiled regex will be faster
          List<String> nmods = i+1<items.size()?Strs.split(items.get(i+1), Strs.PTN_SIGNS, true):null;
          if (backbone.length()==0)
          {
            // N-term mod
            for (String nmod : nmods)
              m.addModificationMatch(ModAttachment.N_TERM, Double.valueOf(nmod));
          }
          else if (i+1<items.size() && nmods!=null)
          {
            // Side-chain mod
            for (String nmod : nmods)
              m.addModificationMatch(backbone.length()-1, Double.valueOf(nmod));
          }
          nmods=(List )Tools.dispose(nmods);
        }
        catch (Exception e)
        {
          e.printStackTrace();
        }
      }
      items=(List )Tools.dispose(items);
      backbone=Tools.dispose(backbone);

//      // deposit the new match for later calls
//      mods_match.put(p, m);

      return m;
    }
    return null;
  }
//  // based on the output in the CSV from mzid file
//  public static PeptideMatch fromMZID(String backbone, String mods)
//  {
//    if (!Strs.isSet(backbone)) return null;
//
//    PeptideMatch m = new PeptideMatch(backbone);
//
//    // Carbamidomethyl:0;iTRAQ4plex114:14;iTRAQ4plex114:15
//    // DVLTLQLEVLMETDSRLHFKIK	iTRAQ4plex114:0;Oxidation:11;iTRAQ4plex114:20;iTRAQ4plex114:22
//    if (Strs.isSet(mods))
//    {
//      String[] items = Strs.split(mods, ';', true);
//      for (String mod : items)
//      {
//        String[] tags = Strs.split(mod, ':');
//        int       pos = Integer.valueOf(tags[1]);
//        Optional<UnimodMod> M = UnimodManager.getModification(tags[0]);
//
//        if (M.isPresent())
//          if (pos==0) m.addModificationMatch(ModAttachment.N_TERM, M.get());
//          else        m.addModificationMatch(pos-1, M.get());
//      }
//    }
//    return m;
//  }

  public static PeptideMatch addScore(PeptideMatch m, Map<String, Object> info, String... keys)
  {
    if (m!=null && Tools.isSet(info) && Tools.isSet(keys))
      for (String key : keys)
        if (info.containsKey(key) && Stats.toDouble(info.get(key).toString())!=null &&
            !Stats.toDouble(info.get(key).toString()).isNaN() && !Stats.toDouble(info.get(key).toString()).isInfinite())
          m.addScore(key, Stats.toDouble(info.get(key).toString()));
    return m;
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
  public static int ppm2bins(double mz0, double mz1, double ppm)
  {
    return (int )(Math.log(1+(mz1-mz0)/mz0)/Math.log(1d+1E-6*ppm));
  }

  /** FragmentAnnotator is a facade that highly simplifies and hides the complexity of annotating a PeakList given a key.
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
    // the fragmenter is responsible for generating a theoretical spectrum from the key
    PeptideFragmenter fragmenter =
        new PeptideFragmenter(EnumSet.of(IonType.b, IonType.y), PeakList.Precision.DOUBLE);

    // the annotator needs to delegate to the fragmenter the fragmentation process and to the internal aligner the tolerance for aligning peaks
    PeptideFragmentAnnotator annotator = new PeptideFragmentAnnotator(fragmenter, new AbsoluteTolerance(0.1));

    // the key as a source for theoretical fragments with annotations
    // Peptide key = Peptide.parse("QVHPDTGISSK");

    // start process: fragmenting key, aligning to the spectrum and reporting matched ions back to the spectrum
    return annotator.annotate(spectrum, peptide, 3);
  }
  public static PeptideSpectrum annotate(PeptideFragmentAnnotator annotator, MsnSpectrum spectrum, PeptideMatch match)
  {
    // convert the spectrum to the suitable format for the annotator
    Peptide      peptide = match.toPeptide(sNumModResolver);
    PeptideSpectrum msms = new PeptideSpectrum(peptide, spectrum.getPrecursor().getCharge());
    for (int i=0; i<spectrum.size(); i++)
      msms.add(spectrum.getMz(i), spectrum.getIntensity(i));

    // start process: fragmenting key, aligning to the spectrum and reporting matched ions back to the spectrum
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
      Multimap<SpectrumIdentifier, PeptideMatch> Bs, int tops, String variety, String score, String delta_score)
  {
    if (!Tools.isSet(Bs)) return As;
    if (As==null) As = HashMultimap.create();

    System.out.print(Bs.size() + " total PSMs from " + Bs.keySet().size() + " MS/MS");
    // combine the PSMs
    DesendScorePeptideMatch sorter = new DesendScorePeptideMatch(score);
    List<PeptideMatch>                   matches = new ArrayList<>();
    Map<String, PeptideMatch>          seq_match = new HashMap<>();
    for (SpectrumIdentifier id : Bs.keySet())
    {
      // System.out.println(" --> " + As.size() + " total PSMs from " + As.keySet().size() + "MS/MS");
      // tag the PSM to the search condition and keeping a map to sort out the PSM by their distinct sequence
      seq_match.clear();
      for (PeptideMatch m : Bs.get(id))
      {
        m.addScore(variety, 1);
        seq_match.put(PSMs.toNumModSequence(m), m);
      }
      // pool the matches for the scan together
      matches.clear(); matches.addAll(seq_match.values());
      if (As.containsKey(id))
        for (PeptideMatch m : As.get(id))
        {
          // making sure this is a new or better hit
          String modseq = PSMs.toNumModSequence(m);
          if (!seq_match.containsKey(modseq) ||
               seq_match.get(modseq).getScore(score)<m.getScore(score)) matches.add(m);
        }
//      // rank the matches by the canonical score in desending order
//      Collections.sort(matches, sorter);
//      // re-assign the rank
//      // TODO need to update the delta score as well
//      for (int i=0; i< matches.size(); i++)
//      {
//        matches.get(i).setRank(i+1);
//        if (Strs.isSet(delta_score) && i<matches.size()-1)
//        {
//          matches.get(i).getScoreMap().remove(delta_score);
//          matches.get(i).addScore(delta_score, matches.get(i).getScore(score)-matches.get(i+1).getScore(score));
//        }
//      }
//      // keep only the top few to preserve the memory footprint
//      As.removeAll(id);
//      As.putAll(id, matches.subList(0, tops <= matches.size() ? tops : matches.size()));

      // keep only the top few to preserve the memory footprint
      As.removeAll(id);
      As.putAll(id, trimByRank(matches, sorter, tops, score, delta_score));
    }
    System.out.println(" --> " + As.size() + " total PSMs from " + As.keySet().size() + " MS/MS");

    return As;
  }
  public static ListMultimap<SpectrumIdentifier, PeptideMatch> add(
      ListMultimap<SpectrumIdentifier, PeptideMatch> id_match, SpectrumIdentifier id, PeptideMatch match,
      int lowest_rank, PSMs.DesendScorePeptideMatch sorter)
  {
    List<PeptideMatch> mm = id_match.get(id);
    if (Tools.isSet(mm) && mm.size()==lowest_rank && sorter.compare(mm.get(mm.size()-1),match)>0)
    {
      if (mm.size()==1) id_match.removeAll(id); else mm.remove(mm.size()-1);
    }
    if (!id_match.containsKey(id) || id_match.get(id).size()<lowest_rank)
    {
      id_match.put(id, match);
      if (id_match.get(id).size()>1)
        Collections.sort(id_match.get(id), sorter);
    }

    return id_match;
  }
  // TODO not WORKING at all. No error emitted but no values got removed either!!!
  public static Multimap<SpectrumIdentifier, PeptideMatch> trimByRank(Multimap<SpectrumIdentifier, PeptideMatch> matches, int lowest_rank)
  {
    if (lowest_rank>0)
      for (SpectrumIdentifier id : matches.keySet())
      {
        Collection<PeptideMatch> mm = matches.get(id);
        Iterator<PeptideMatch> itr = mm.iterator();
        while (itr.hasNext())
        {
          PeptideMatch m = itr.next();
          if (m.getRank()>lowest_rank) itr.remove();
        }
      }

    return matches;
  }
  public static Collection<PeptideMatch> trimByRank(
      List<PeptideMatch> matches, DesendScorePeptideMatch sorter, int tops, String score, String delta_score)
  {
    if (tops==0) return matches;

    // rank the matches by the canonical score in desending order
    Collections.sort(matches, sorter);
    // re-assign the rank
    // TODO need to update the delta score as well
    for (int i=0; i< matches.size(); i++)
    {
      matches.get(i).setRank(i+1);
      if (Strs.isSet(delta_score) && i<matches.size()-1)
      {
        matches.get(i).getScoreMap().remove(delta_score);
        matches.get(i).addScore(delta_score, matches.get(i).getScore(score)-matches.get(i+1).getScore(score));
      }
    }
    return matches.subList(0, tops <= matches.size() ? tops : matches.size());
  }
  public static PeptideMatch clone(PeptideMatch m)
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
  public static PeptideMatch clone2LOD(PeptideMatch m, double multiple)
  {
    // look for the lowest rank for this engine
    PeptideMatch dummy = PSMs.clone(m); dummy.getScoreMap().clear();
    // set the score multiples
    for (String score : m.getScoreMap().keySet())
      dummy.addScore(score, m.getScore(score)*multiple);

    return dummy;
  }
  public static PeptideMatch worst(Engine E, Collection<PeptideMatch> matches)
  {
    PeptideMatch worst=null;
    for (PeptideMatch match : matches)
      if (worst==null || (E.compareTo(worst.getScoreMap(), match.getScoreMap())<0)) worst=match;

    return worst;
  }

  public static Multimap<String, PeptideMatch> toBackboneMatchMap(Collection<PeptideMatch> matches)
  {
    if (matches==null) return null;

    Multimap<String, PeptideMatch> seq_match = HashMultimap.create();
    for (PeptideMatch match : matches)
      seq_match.put(match.toSymbolString(), match);

    return seq_match;
  }
  // Convert a q-val or other probablistic score to a canonical score where higher is better in linear scale
  public static List<PeptideMatch> setCanonicalDeltaScoresByQ(Collection<PeptideMatch> hits, String qval)
  {
    if (!Tools.isSet(hits)) return null;

    List<PeptideMatch> psm = new ArrayList<>(hits);
    // look for the smallest value above the zero
    double floor=Double.MAX_VALUE;
    for (PeptideMatch m : psm)
      if (m.getScore(qval)>0 && m.getScore(qval)<floor) floor=m.getScore(qval);

    for (PeptideMatch m : psm)
    {
      m.addScore(SCR_CANONICAL, -10d*Math.log10(m.getScore(qval)>0?m.getScore(qval):floor/2d));
    }
    Collections.sort(psm, new DesendScorePeptideMatch(SCR_CANONICAL));
    // set the delta scores
    if (psm.size()>1)
      for (int i=0; i<psm.size()-1; i++)
      {
        psm.get(i).addScore(SCR_DELTA, psm.get(i).getScore(SCR_CANONICAL)-psm.get(i+1).getScore(SCR_CANONICAL));
      }
    psm.get(psm.size()-1).addScore(SCR_DELTA, 0);


    return psm;
  }
/*
  */
/**@param mm are the input key matches
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
