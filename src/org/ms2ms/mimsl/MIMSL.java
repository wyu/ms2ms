package org.ms2ms.mimsl;

import com.google.common.collect.Range;
import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.expasy.mzjava.proteomics.mol.AAMassCalculator;
import org.expasy.mzjava.proteomics.ms.consensus.PeptideConsensusSpectrum;
import org.expasy.mzjava.proteomics.ms.spectrum.PepFragAnnotation;
import org.expasy.mzjava.proteomics.ms.spectrum.PepLibPeakAnnotation;
import org.ms2ms.algo.Peaks;
import org.ms2ms.mzjava.AnnotatedPeak;
import org.ms2ms.mzjava.AnnotatedSpectrum;
import org.ms2ms.nosql.ms.HBaseProteomics;
import org.ms2ms.algo.MsStats;
import org.ms2ms.utils.Tools;

import java.io.IOException;
import java.util.*;

/** The core algorithms for Mrm inspired MS/MS spectral lookup
 *
 * Created by wyu on 4/22/14.
 */
public class MIMSL
{
  static class ScoreDesendComparator implements Comparator<AnnotatedSpectrum>
  { public int compare(AnnotatedSpectrum o1, AnnotatedSpectrum o2) { return o1!=null && o2!=null ?
    Double.compare(o2.getScore(AnnotatedSpectrum.SCR_MIMSL), o1.getScore(AnnotatedSpectrum.SCR_MIMSL)):0; } }

  public static void randHBasePeakList(Random rand)
  {

  }
  synchronized public static List<AnnotatedSpectrum> run(Peak[] precursors, MimslSettings settings, Peak... frags) throws IOException
  {
    long nsec = System.nanoTime();

    List<AnnotatedSpectrum> candidates = new ArrayList<AnnotatedSpectrum>();
    candidates.addAll(setStatus(HBaseProteomics.query(precursors, settings, 0d, frags), PeptideConsensusSpectrum.Status.NORMAL));
    // add 7 da offset to simulate decoy matches since this is not a common offset due to mod or mutation
    candidates.addAll(setStatus(HBaseProteomics.query(precursors, settings, 7d, frags), PeptideConsensusSpectrum.Status.DECOY));
    System.out.println("Query time: " + Tools.d2s(1E-9*(System.nanoTime()-nsec), 2) + " sec.");

    // calculate the score by hypergeometric model
    candidates = (List<AnnotatedSpectrum> )score(candidates, settings.getFragmentTol());

    nsec = System.nanoTime();
    HBaseProteomics.loadPeakLists(candidates);

    fdr(candidates);
    System.out.println("Load time: " + Tools.d2s(1E-9*(System.nanoTime()-nsec), 2) + " sec.");

    return candidates;
  }
  @Deprecated
  synchronized public static List<AnnotatedSpectrum> run(PeakList<PepLibPeakAnnotation> ions, byte[] spec_type, Tolerance precursor, Tolerance frag) throws IOException
  {
    long nsec = System.nanoTime();

    List<AnnotatedSpectrum> candidates = new ArrayList<AnnotatedSpectrum>();
    candidates.addAll(setStatus(HBaseProteomics.query(ions, spec_type, precursor, 0d), PeptideConsensusSpectrum.Status.NORMAL));
    // add 7 da offset to simulate decoy matches since this is not a common offset due to mod or mutation
    candidates.addAll(setStatus(HBaseProteomics.query(ions, spec_type, precursor, 7d), PeptideConsensusSpectrum.Status.DECOY));
    System.out.println("Query time: " + Tools.d2s(1E-9*(System.nanoTime()-nsec), 2) + " sec. m/z" +
        Tools.d2s(precursor.getMin(ions.getPrecursor().getMz()),4) + " to " + Tools.d2s(precursor.getMax(ions.getPrecursor().getMz()),4));

    // calculate the score by hypergeometric model
    candidates = (List<AnnotatedSpectrum> )score(candidates, frag);

    nsec = System.nanoTime();
    HBaseProteomics.loadPeakLists(candidates);

    fdr(candidates);
    System.out.println("Load time: " + Tools.d2s(1E-9*(System.nanoTime()-nsec), 2) + " sec.");

    return candidates;
  }
  @Deprecated
  public static PepLibPeakAnnotation getSignatureAnnotation(double mz, Collection<PepLibPeakAnnotation> annos, double min_mz, double precursor_mz)
  {
    if (mz>0 && annos!=null)
      for (PepLibPeakAnnotation anno : annos)
      {
        // TODO need to store the calculated m/z of the fragment with anno
        // some of the published transition use 2+ ions below precursor m/z
        PepFragAnnotation frag = anno.getOptFragmentAnnotation().get();
        // for a default value of 1+ if the charge is not set
        boolean OK = AAMassCalculator.getInstance().calculateNeutralMolecularMass(mz, frag.getCharge()==0?1:frag.getCharge()) > precursor_mz &&
          mz > min_mz && Peaks.isType(anno, IonType.b, IonType.y) && !Peaks.isType(anno, IonType.p, IonType.unknown) && frag.getNeutralLoss().getMolecularMass()==0d;
        if (OK) return anno;
      }

    return null;
  }
  public static PepLibPeakAnnotation getSignatureAnnotation(Collection<PepLibPeakAnnotation> annos)
  {
    if (annos!=null)
      for (PepLibPeakAnnotation anno : annos)
        if (Peaks.isType(anno, IonType.b, IonType.y) && !Peaks.isType(anno, IonType.p, IonType.unknown) &&
            anno.getOptFragmentAnnotation().get().getNeutralLoss().getMolecularMass()==0d &&
          anno.getOptFragmentAnnotation().get().getIsotopeCount()==0) return anno;

    return null;
  }

  /** Extract the signature peaks from an annotated MS/MS.
   *
   *  !! invalidate the noise peaks a prior by min_merge_count or local noise model !!
   *
   * @param msms
   * @param half_width
   * @return
   */
  public static List<AnnotatedPeak> toSignature(PeakList msms, double half_width, double min_mz, int tops, double min_snr)
  {
    if (msms == null || msms.size()==0 || msms.getAnnotationIndexes().length==0) return null;
    // setup a map to avoid duplicated signature of the same m/z
    SortedMap<Double, AnnotatedPeak> mz_signature = new TreeMap<Double, AnnotatedPeak>();
    // working objects
    List<AnnotatedPeak> orphans = new ArrayList<AnnotatedPeak>();
    AnnotatedPeak         below = null; // the lengest frag just below the precursor m/z

    // only step thro the indices with peak annotation
    for (int i :msms.getAnnotationIndexes())
    {
      PepLibPeakAnnotation OK = getSignatureAnnotation(msms.getAnnotations(i));
      if (OK!=null)
      {
        PepFragAnnotation f = OK.getOptFragmentAnnotation().get();
        double    pk_counts = Peaks.countValid( msms, msms.getMz(i)-half_width, msms.getMz(i)+half_width),
                       base = Peaks.getBaseline(msms, msms.getMz(i)-half_width, msms.getMz(i)+half_width, tops, true),
                         mh = AAMassCalculator.getInstance().calculateNeutralMolecularMass(msms.getMz(i), f.getCharge()==0?1:f.getCharge());
        AnnotatedPeak   ion = new AnnotatedPeak(f.getTheoreticalMz(), msms.getIntensity(i), f.getCharge()==0?1:f.getCharge(), base!=0?Math.abs(msms.getIntensity(i)/base):-1d);

        // make sure the frag is more intense than the baseline
        if (mh>=min_mz && (msms.getPrecursor().getCharge()==1 || mh>msms.getPrecursor().getMz()))
        {
          if      (pk_counts>0 && ion.getSNR()>=min_snr)
          {
            if (mz_signature.get(ion.getMz())==null ||
                mz_signature.get(ion.getMz()).getIntensity()<ion.getIntensity()) mz_signature.put(ion.getMz(),ion);
          }
          else if (pk_counts<5) orphans.add(ion);
        }
        if (msms.getPrecursor().getCharge()>1 &&
            msms.getMz(i)<msms.getPrecursor().getMz()-28d &&
            ion.getSNR()>=min_snr && (below==null || below.getMz()<ion.getMz())) below = ion;
      }
    }

    if (mz_signature.size() < tops && below!=null) mz_signature.put(below.getMz(), below);
    if (mz_signature.size() < tops && Tools.isSet(orphans))
      for (AnnotatedPeak p : orphans)
        mz_signature.put(p.getMz(), p);

    return new ArrayList<AnnotatedPeak>(mz_signature.values());
  }
  public static Range<Double> range(double mz, Tolerance tol)
  {
    return Range.closed(tol.getMin(mz), tol.getMax(mz));
  }
  public static Collection<Range<Peak>> range(Peak p, Tolerance tol, int... zs)
  {
    Collection<Range<Peak>> ranges = new ArrayList<Range<Peak>>();
    // if alternative not present, use the main peak
    if (!Tools.isSet(zs)) zs = new int[] {p.getCharge()};
    // step thro the charge states
    for (int z : zs)
    {
      double mz = 1.00078 + Peaks.toMass(p) / (double) z;
      Peak p1 = new Peak(p), p2 = new Peak(p);
      p1.setMzAndCharge(tol.getMin(mz), z);
      p2.setMzAndCharge(tol.getMax(mz), z);
      ranges.add(Range.closed(p1, p2));
    }
    return ranges;
  }
  public static Range<Double> range(double mass, int z, Tolerance tol)
  {
    double mz = (1.00078 + mass) / (double )z;
    return Range.closed(tol.getMin(mz), tol.getMax(mz));
  }
  // Given a list of precursor mz's, enumerate the m/z slices that are consistent with the mz tolerance and isolation width
  public static Collection<Range<Peak>> enumeratePrecursors(Tolerance tol, int zfloat, Peak... precursors)
  {
    if (!Tools.isSet(precursors)) return null;

    List<Range<Peak>> slices = new ArrayList<Range<Peak>>();
    // expand the precursors into m/z slices
    for (Peak ion : precursors)
    {
      int zs[] = null;
      if (ion.getCharge()!=0 && zfloat>0)
      {
        zs = new int[zfloat*2+1];
        zs[0] = ion.getCharge();
        // check the neighboring charges if asked
        for (int z = 1; z < zfloat; z++)
        {
          zs[(z-1)*2+1]=ion.getCharge()-z;
          zs[(z-1)*2+2]=ion.getCharge()+z;
        }
      }
      slices.addAll(range(ion, tol, zs));
    }
    return Tools.merge(slices);
  }
  public static <T extends PeptideConsensusSpectrum> Collection<T> setStatus(Collection<T> spectra, PeptideConsensusSpectrum.Status status)
  {
    if (Tools.isSet(spectra))
      for (T spec : spectra) spec.setStatus(status);

    return spectra;
  }
  public static Collection<AnnotatedSpectrum> score(Collection<AnnotatedSpectrum> candidates, Tolerance tol)
  {
    if (!Tools.isSet(candidates)) return candidates;

    for (AnnotatedSpectrum spec : candidates)
    {
      long bins = Math.round(1000 / (tol.getMax(spec.getPrecursor().getMz())- tol.getMin(spec.getPrecursor().getMz())));
      int nmatch = spec.getIonMatched(), nfrag = spec.getIonQueried(), nsig = spec.getIonIndexed();
      // long success, long trials, long success_population, long population
      spec.setScore(AnnotatedSpectrum.SCR_MIMSL,       -1d * MsStats.hypergeometricPval1(nmatch, nfrag, nsig, bins));
      if (nsig > 1) nsig--; else if (nfrag > 1) nfrag--;
      spec.setScore(AnnotatedSpectrum.SCR_MIMSL_DELTA, -1d * MsStats.hypergeometricPval1(nmatch, nfrag, nsig, bins) - spec.getScore(AnnotatedSpectrum.SCR_MIMSL));
    }
    return candidates;
  }
  public static List<AnnotatedSpectrum> fdr(List<AnnotatedSpectrum> candidates)
  {
    if (Tools.isSet(candidates))
    {
      // rank the candidates by their scores
      Collections.sort(candidates, new ScoreDesendComparator());
      // TODO fdr calculation
    }
    return candidates;
  }
  public static StringBuffer printCandidates(StringBuffer buf, Collection<AnnotatedSpectrum> candidates)
  {
    if (!Tools.isSet(candidates)) return buf;
    if (buf==null) buf=new StringBuffer();

    buf.append("score\tdelta\tvotes\tverdict\tdecoy\tppm\tPeptide\tm/z\tz\tSig\tunmatch\tprotein\n");

    for (AnnotatedSpectrum candidate : candidates)
    {
      String[] peptide = candidate.getComment().split("\\^");
      boolean decoy = (candidate.getStatus().equals(PeptideConsensusSpectrum.Status.DECOY));
      buf.append(Tools.d2s(candidate.getScore(AnnotatedSpectrum.SCR_MIMSL),       2) + "\t");
      buf.append(Tools.d2s(candidate.getScore(AnnotatedSpectrum.SCR_MIMSL_DELTA), 2) + "\t");
      buf.append(candidate.getIonMatched() + "\t");
      buf.append(candidate.getStatus() + "\t");
      buf.append((decoy?"---":Tools.d2s(Peaks.toPPM(candidate.getPrecursor().getMz(), candidate.getMzQueried()), 2)) + "\t");
      buf.append((decoy?"---":peptide[0]) + "\t");
      buf.append(Tools.d2s(candidate.getPrecursor().getMz(), 4) + "\t");
      buf.append(candidate.getPrecursor().getCharge() + "\t");
      buf.append(candidate.getIonIndexed() + "\t");
      buf.append((candidate.getIonIndexed()-candidate.getIonMatched()) + "\t");
      buf.append((decoy?"---":peptide.length>1?peptide[1]:"") + "\n");
    }
    return buf;
  }

/*  protected static final String GRP_SETTINGS     = "Parameters";
  protected static final String GRP_IONS         = "Signature Ions";
  protected static final String GRP_CANDIDATE    = "";

  public static float        sDecoyOffset = 7f;
  public static boolean      sVerbose = true;
  public static StringBuffer sMessage = new StringBuffer();

  private boolean                 mIsDirty = false;
  private MsSettings mSettings = new MsSettings();
  private Multimap<Integer, MsIon> mFrags = TreeMultimap.create(),
    mPrecursorMzs = TreeMultimap.create();

  protected Float mElipseMsec = 0f;
  private long mBinSize = 0L;
  private int mRowLimit = 10;

  private Analysis mMatches = new Analysis(0f), mDecoys = new Analysis(sDecoyOffset);
  private String mSettingsName = "LTQ 0.75Da^0.35Da", mLibName = MsMsDictionary.LIB_HUMAN;
  private Collection<PTM.ePtmClassification> mPtmInclusions = null;
  private Boolean mPtmArtifects = false, mPtmTr = false, mPtmGlyco = false, mPtmAA = false, mPtmIsotope = false;

  protected MultiTreeMap<Float, PTM> mPTMs;

  public MimpCore()                    { super(); }
  public MimpCore(MsSettings settings) { super(); mSettings = settings; }
  public MimpCore(String lib) { super(); init(lib); }
  public MimpCore(String lib, MsSettings settings) { super(); mSettings = settings;  init(lib); }

  // the Getters
  public Analysis getMatches() { return mMatches; }
  public Analysis getDecoys() { return mDecoys; }

  public boolean isDirty() { return mIsDirty; }
  public MsSettings getSettings()    { return mSettings; }
  public    Multimap<Integer, MsIon> getPrecursorMzs() { return mPrecursorMzs; }
  public    Multimap<Integer, MsIon> getFragmentMzs() { return mFrags; }
  //public SortedMap<PeptideHit, MsMsAssignment>    getCandidates()  { return mCandidates; }
  //protected MultiTreeMapSet<Long, Float> getMatches()     { return mMatches; }
  protected MultiTreeMap<Float, PTM> getPtmpMap() { return mPTMs; }
  protected Collection<Float> getPtmMs() { return mPTMs != null ? mPTMs.keySet() : null; }

  public int getRowLimit() { return mRowLimit; }
  public MimpCore setRowLimit(int s) { mRowLimit = s; return this; }

  public MimpCore isDirty(boolean s) { mIsDirty = s; return this; }
  public MimpCore setPrecursor(Double s, Integer z)
  {
    if (mPrecursorMzs != null) mPrecursorMzs.clear(); else mPrecursorMzs = TreeMultimap.create();
    mPrecursorMzs.put(z, new MsPeak(s, 0d, z)); return this;
  }
  public MimpCore addPrecursor(MsIon ion)
  {
    if (mPrecursorMzs == null) mPrecursorMzs = TreeMultimap.create();
    // make sure the precursor is not already deposited
    if (!MsIon_Util.contains(mPrecursorMzs.cells(ion.getCharge()), ion,
      getSettings().getEffectivePrecursorTol(ion.getMz()))) mPrecursorMzs.put(ion.getCharge(), ion);

    return this;
  }
  public MimpCore setPrecursors(Double... mz)
  {
    if (mPrecursorMzs != null) mPrecursorMzs.clear(); else mPrecursorMzs = TreeMultimap.create();
    for (Double m : mz) mPrecursorMzs.put(0, new MsPeak(m, 0d, 0));
    return this;
  }
  public MimpCore setFrags(Double... mz)
  {
    mFrags.clear();
    for (Double m : mz) mFrags.put(0, new MsPeak(m, 0d, 0));
    return this;
  }
  public MimpCore setFrags(Collection<MsIon> mzs)
  {
    mFrags.clear();
    for (MsIon m : mzs) if (m.isValid()) mFrags.put(m.getCharge(), m);
    return this;
  }
  public Map<PeptideHit, MsMsAssignment> fdr(int row_limit, double fold_x_background)
  {
    if (mMatches.iCandidates == null) throw new RuntimeException("NULL candidate pool!");

    SortedSet<SimpleCandidate> hits = new TreeSet<SimpleCandidate>(SimpleCandidate.SORT_SCORE_KEY_DESEND);
    Iterator<SimpleCandidate> itr = mMatches.iCandidates.iterator();
    while (itr.hasNext())
      itr.next().setVerdict(Candidate.eVerdict.normal);

    hits.addAll( mMatches.iCandidates);

    //hits.addAll(mMatches.iCandidates);
    if (mDecoys.iCandidates != null)
    {
      itr = mDecoys.iCandidates.iterator();
      while (itr.hasNext())
        itr.next().setVerdict(Candidate.eVerdict.decoy);

      hits.addAll(mDecoys.iCandidates);
    }
    if (!Toolbox.isSet(hits)) return null;

    // figure out the score of the best decoy
    double score_decoy = 0, fdr = 200d / (double )hits.size(), ln_folds = Math.log(fold_x_background);
    for (SimpleCandidate H : hits)
      // the hits are sorted by the score in decending order
      if (Toolbox.equals(H.getVerdict(), Candidate.eVerdict.decoy))
      //{ score_decoy = H.getScore() + H.getScoreDelta(); break; }
      { score_decoy = H.getScore() + ln_folds; break; }

    Map<PeptideHit, MsMsAssignment> candidates = new TreeMap<PeptideHit, MsMsAssignment>(PeptideHit.SORT_SCORE_DEVI);
    for (SimpleCandidate H : hits)
    {
      candidates = Toolbox.equals(H.getVerdict(), Candidate.eVerdict.decoy) ?
        mDecoys.score(candidates, H.getPrimaryKey(), true) :
        mMatches.score(candidates, H.getPrimaryKey(), false);

      if (!Toolbox.equals(H.getVerdict(), Candidate.eVerdict.accepted) && candidates.size() >= row_limit) break;
    }

    int accepted = 0;
    for (PeptideHit H : candidates.keySet())
    {
      H.setScoreDelta(null);
      H.removeProperty("ScoreDelta");
      H.removeProperty("dS");

      // the hits are sorted by the score in decending order
      if (H.getScore() >= score_decoy)
      {
        H.setVerdict(Candidate.eVerdict.accepted);
        accepted++;
      }
      else
        H.setVerdict(Candidate.eVerdict.mistake);
    }
    List<PeptideHit> candidates_ordered = new ArrayList<PeptideHit>(candidates.keySet());
    PeptideCandidate_Util.calcScoreDelta(candidates_ordered, null, false); // all charges

    yell("score threshold: " + Toolbox.d2s(score_decoy, 3) + " @ fdr of " + Toolbox.d2s(fdr, 3) + "%");
    yell("counts (accepted/total): " + accepted + "/" + hits.size());

    return candidates;
  }
  public boolean isReady()
  {
    return Toolbox.isSet(mFrags) && (MsMsDictionary.isReady(mLibName));
  }
  public boolean canProceed()
  {
    if (!Toolbox.isSet(mPrecursorMzs) && (!Toolbox.isSet(mFrags) || mFrags.size() < 3))
      throw new UserException("Must have at least 3 fragments if the precursor is absent!");

    if (Toolbox.isSet(mPTMs))
    {
      if (!Toolbox.isSet(mPrecursorMzs))
        throw new UserException("Precursor required for Mod-tolerant Query!");
      if (getSettings().getEffectivePrecursorTol(500d) > 0.1d)
        throw new UserException("High-resolution Precursor required for Mod-tolerant Query!");

      for (MsIon ion : mPrecursorMzs.values())
        if (ion.getCharge() == null || ion.getCharge() == 0)
          throw new UserException("Precursor Charge required for Mod-tolerant Query!");

      if (!Toolbox.isSet(mFrags) || mFrags.size() < 3)
        throw new UserException("Must have at least 3 fragments for Mod-tolerent Query!");
    }
    return true;
  }
  private MultiTreeMap<Float, PTM> populatePTMs()
  {
    mPtmInclusions = new HashSet<PTM.ePtmClassification>();
    if (mPtmArtifects) mPtmInclusions.add(PTM.ePtmClassification.artifact);
    if (mPtmAA)
    {
      mPtmInclusions.add(PTM.ePtmClassification.aa);
      mPtmInclusions.add(PTM.ePtmClassification.residue);
    }
    if (mPtmGlyco)
    {
      mPtmInclusions.add(PTM.ePtmClassification.glyco);
      mPtmInclusions.add(PTM.ePtmClassification.nlink);
      mPtmInclusions.add(PTM.ePtmClassification.olink);
    }
    if (mPtmTr)
    {
      mPtmInclusions.add(PTM.ePtmClassification.pre_tr);
      mPtmInclusions.add(PTM.ePtmClassification.co_tr);
      mPtmInclusions.add(PTM.ePtmClassification.post_tr);
    }
    if (mPtmIsotope)
    {
      mPtmInclusions.add(PTM.ePtmClassification.chemical);
      mPtmInclusions.add(PTM.ePtmClassification.isotope);
      mPtmInclusions.add(PTM.ePtmClassification.peptide_syn);
    }
    mPTMs = Toolbox.isSet(mPtmInclusions) ? Sequence_Util.newUnimodPTMs(mPtmInclusions) : null;
    return mPTMs;
  }
  synchronized public boolean run() throws Exception
  {
    yell("{code}");

    mElipseMsec = (float )System.nanoTime();
    //mBinSize    = (long )(1000 / mSettings.getEffectivePrecursorTol(
    //    Toolbox.isSet(mPrecursorMzs) ? Toolbox.front(mPrecursorMzs.values()).getMz() : 600f));
    mBinSize    = (long )(1000 / mSettings.getEffectiveFragmentTol(
      Toolbox.isSet(mFrags) ? MsIon_Util.meanX(mFrags.values()) : 600f));

    mMatches.run();
    mDecoys.run();
    MsMsDictionary.cells(mLibName).close(); isDirty(false);

    yell("{code}");

    //dictionary.announce(mimp.getCandidates(), "The final candidates", 10);
    return true;
  }
  protected void init(String lib)
  {
    mLibName  = lib;
    mSettings = MsSettings.LTQ;
    // null PTMs for now
    MsMsDictionary.prepare(mLibName, false, null); // not wait for the finish
  }

  // Given a list of precursor mz's and potential PTM, enumerate the m/z slices that are consistent with the mz tolerance and isolation width
  public static Collection<Range<Float>> enumeratePrecursor(MsSettings settings, Collection<Float> mods_, Collection<MsIon> precursors, Integer z)
  {
    if (precursors == null || settings == null || z == null || z == 0) return null;

    List<Range<Float>> slices = new ArrayList<Range<Float>>();

    // need to have at least one member
    Collection<Float> mods = new ArrayList<Float>();  mods.add(0f);
    if (mods_ != null) mods.addAll(mods_);

    // make sure the isolation window is set
    Range<Float> isolation = settings.getPrecursorMzRange() != null ? settings.getPrecursorMzRange() : null;

    // expand the precursors into m/z slices
    for (MsIon ion : precursors)
    {
      float delta = (float )settings.getEffectivePrecursorTol(ion.getMz());
      for (Float offset : mods)
      {
        float m = (float )(ion.getMz() - (offset / z));
        if (isolation != null)
        {
          float mid = 0f;
          while (mid <= 0) // always true
          {
            if (mid + delta <= isolation.getUpper())
              slices.add(new Range<Float>(m + Math.max(mid - delta, isolation.getLower()),
                m + Math.min(mid + delta, isolation.getUpper())));
            mid -= (Mass.HYDROGEN_MONO / z);
            if (mid - delta < isolation.getLower()) break;
          }

          mid = (float )(Mass.HYDROGEN_MONO / z);
          while (mid >= 0)
          {
            if (mid - delta >= isolation.getLower())
              slices.add(new Range<Float>(m + Math.max(mid - delta, isolation.getLower()),
                m + Math.min(mid + delta, isolation.getUpper())));
            mid += (Mass.HYDROGEN_MONO / z);
            if (mid + delta > isolation.getUpper()) break;
          }
        }
        else
        {
          slices.add(new Range<Float>(m - delta, m + delta));
        }
      }
    }
    // merge the slices
    Collections.sort(slices);
    for (int i = 0; i < slices.size(); i++)
      if (slices.cells(i).isValid())
        for (int j = i+1; j < slices.size(); j++)
          if (slices.cells(j).isValid() && slices.cells(i).hasOverlap(slices.cells(j)))
          {
            slices.cells(i).extend(slices.cells(j)); slices.cells(j).invalidate();
          }

    Iterator<Range<Float>> itr = slices.iterator();
    while (itr.hasNext()) if (!itr.next().isValid()) itr.remove();

    return slices;
  }
  // making the call whether the candidates qualified as a match
  public boolean isMatch(Map<PeptideHit, MsMsAssignment> candidates,
    Double min_score, Integer min_rank, Double min_delta)
  {
    if (min_rank == null) min_rank = 1;
    if (!Toolbox.isSet(candidates)) return false;

    int order = 0;
    for (PeptideHit H : candidates.keySet())
    {
      XYPoint pt = calcResidual(H.getMassDeviation().floatValue(), H.getCharge());
      if (H.is(Candidate.eVerdict.accepted) &&
        pt.getX() == 0 &&
        (min_delta == null || H.getScoreDelta() == null || H.getScoreDelta() > min_delta) &&
        (min_score == null || H.getScore() >= min_score))
      {
        return true;
      }
      if (++order > min_rank) break;
    }
    return false;
  }

  public Map<PeptideHit, MsMsAssignment> announce(String title, int rows, double folds_x_background)
  {
    return announce(title, fdr(rows, folds_x_background));
  }
  public Map<PeptideHit, MsMsAssignment> announce(String title, Map<PeptideHit, MsMsAssignment> candidates)
  {
    if (!Toolbox.isSet(candidates)) return null;

    System.out.println("\n" + title);
    System.out.println("||score||delta||votes||verdict||decoy||offset||ppm||()||Peptide||()||Charge||m/z||ID||Sig||Unmatched||");

    // (K)AAC[Carbamidomethylation]LLPKLDELR(D), +3, m/z466.931
    int order = 0;
    for (PeptideHit H : candidates.keySet())
    {
      XYPoint pt = calcResidual(H.getMassDeviation().floatValue(), H.getCharge());
      System.out.println("|" + Toolbox.d2s(H.getScore(), 1) + "|" + Toolbox.d2s(H.getScoreDelta(), 2) + " |" +
        H.getRank() + "|" + H.getVerdict() + "|" + (H.isDecoy() ? "Y" : " ") + "|" +
        Toolbox.d2s(pt.getX(), 0) + "|" + Toolbox.d2s(1E6 * pt.getY() / H.getMzObs(), 1) + "|" +
        H.getPrevAA() + "|" + Wiki_Util.getBackboneWiki(H) + "|" + H.getNextAA() + "|+" + H.getCharge() + "|" +
        Toolbox.d2s(H.getMzObs(), 3) + "|" + H.getPeptideHitId().toString() + "|" +
        Toolbox.d2s(H.getDetectibility(), 0) + "|" +  H.getDeNovoCandidate() + "|");
      //if (++order > rows) break;
    }
    return candidates;
  }
  public List<PeptideHit> announce(String title, Map<PeptideHit, MsMsAssignment> candidates, int row_limit, boolean stop_at_decoy)
  {
    List<PeptideHit> output = new ArrayList<PeptideHit>();

    if (!Toolbox.isSet(candidates)) return output;

    System.out.println("\n" + title);
    System.out.println("||score||delta||votes||verdict||decoy||offset||ppm||()||Peptide||()||Charge||m/z||ID||Sig||Unmatched||");

    int order = 0;
    for (PeptideHit H : candidates.keySet())
    {
      XYPoint pt = calcResidual(H.getMassDeviation().floatValue(), H.getCharge());
      System.out.println("|" + Toolbox.d2s(H.getScore(), 1) + "|" + Toolbox.d2s(H.getScoreDelta(), 2) + " |" +
        H.getRank() + "|" + H.getVerdict() + "|" + (H.isDecoy() ? "Y" : " ") + "|" +
        Toolbox.d2s(pt.getX(), 0) + "|" + Toolbox.d2s(1E6 * pt.getY() / H.getMzObs(), 1) + "|" +
        H.getPrevAA() + "|" + Wiki_Util.getBackboneWiki(H) + "|" + H.getNextAA() + "|+" + H.getCharge() + "|" +
        Toolbox.d2s(H.getMzObs(), 3) + "|" + H.getPeptideHitId().toString() + "|" +
        Toolbox.d2s(H.getDetectibility(), 0) + "|" +  H.getDeNovoCandidate() + "|");

      // append the title
      H.setDescription(Toolbox.extend("{" + title + "}", H.getDescription(), ", "));
      output.add(H);
      // stop if this is a decoy hit,
      if (output.size() > row_limit || (stop_at_decoy && Toolbox.equals(H.getVerdict(), Candidate.eVerdict.decoy))) break;
    }
    return output;
  }
  // prepare for the HTML output to be displayed in a web page
  public HTMLTag toHTML(String title, Map<PeptideHit, MsMsAssignment> candidates, boolean brief)
  {
    Div report = new Div().setId("id_report"), headline = report.addDiv().setClass("cls_title");

    if (!Toolbox.isSet(candidates))
    {
      headline.addDiv(new Div("No matching candidate can be found."));
      headline.addDiv(new Div("Candidate: "  + (mMatches.iCandidates       == null ? " NULL" : mMatches.iCandidates.size() + "")));
      headline.addDiv(new Div("Precursors: " + (mMatches.iPrecursorMatches == null ? " NULL" : mMatches.iPrecursorMatches.size() + "")));

      headline.addDiv(new Div("Library name: "  + getLibName()));
      headline.addDiv(new Div("Libraries: "  + MsMsDictionary.sLibrary.size() + ", " + Toolbox.toString(MsMsDictionary.sLibrary.keySet(), ",")));

      MsMsDictionary dic = MsMsDictionary.sLibrary.cells(getLibName());
      if (dic != null)
      {
        headline.addDiv(new Div("Precursor Index: "  + (dic.mPrecursors != null ? dic.mPrecursors.size() : null)));
        headline.addDiv(new Div(" Fragment Index: "  + (dic.mFragments != null ? dic.mFragments.size() : null)));
      }
      else
      {
        headline.addDiv(new Div("NULL dictionary"));
      }

      return report;
    }

    if (Toolbox.isSet(title)) headline.addSpan(title + ", ");

    String url = "/Silo?action=revise", line = null;
    if (Toolbox.isSet(getPrecursorMzs()))
    {
      url += "&precursors="; line = "precursors=";
      for (MsIon ion : getPrecursorMzs().values())
      {
        String item = ion.getMz() + (ion.getCharge() > 0 ? "/" + ion.getCharge() : "");
        url  = Toolbox.extend(url,  item, " ");
        line = Toolbox.extend(line, item, " ");
      }
    }
    if (Toolbox.isSet(getFragmentMzs()))
    {
      url += "&fragments="; line += ", fragments=";
      for (MsIon ion : getFragmentMzs().values())
      {
        String item = ion.getMz() + (ion.getCharge() > 0 ? "/" + ion.getCharge() : "");
        url  = Toolbox.extend(url,  item, " ");
        line = Toolbox.extend(line, item, " ");
      }
    }
    url += "&settings=" + getSettingsStr();

    //headline.addSpan(new Link(url, "revise"));

    Table tbl = report.addTable().setId("id_report_table");

    Tr header = tbl.addRow().setId("id_report_table_header").setClass("cls_row_head");
    if (brief)
    {
      header.addCell("score");
      header.addCell("offset");
      header.addCell("ppm");
      header.addCell("z");
      header.addCell("m/z");
      header.addCell("peptide");
    }
    else
    {
      header.addCell("score");
      header.addCell("delta");
      header.addCell("votes");
      header.addCell("verdict");
      header.addCell("decoy");
      header.addCell("offset");
      header.addCell("ppm");
      header.addCell("z");
      header.addCell("m/z");
      header.addCell("()");
      header.addCell("peptide");
      header.addCell("()");
      header.addCell("ID");
      header.addCell("sig");
      header.addCell("unmatched");
    }

    // (K)AAC[Carbamidomethylation]LLPKLDELR(D), +3, m/z466.931
    int order = 0;
    for (PeptideHit H : candidates.keySet())
    {
      Tr row = tbl.addRow().setClass(++order % 2 == 0 ? "cls_row_shade" : "cls_row_light");

      XYPoint pt = calcResidual(H.getMassDeviation().floatValue(), H.getCharge());
      String peptide = H.getBackboneCutPattern(false).replaceAll("C\\[Carbamidomethylation\\]", "c"
      ).replaceAll("\\[Carbamidomethylation\\]C", "c"
      ).replaceAll("\\[Carbamidomethyl\\]C", "c"
      ).replaceAll("C\\[Carbamidomethyl\\]", "c");
      if (H.isDecoy()) peptide = "<font color='gray'>" + peptide + "</font>";

      if (brief)
      {
        row.addCell(Toolbox.d2s(H.getScore(), 1));
        row.addCell(Toolbox.d2s(pt.getX(), 0));
        row.addCell(Toolbox.d2s(1E6 * pt.getY() / H.getMzObs(), 1));
        row.addCell(H.getCharge() + "");
        row.addCell(Toolbox.d2s(H.getMzObs(), 3));
        Td C = new Td();
        C.addContentWithoutEscaping(peptide);
        row.addCell(C);
      }
      else
      {
        row.addCell(Toolbox.d2s(H.getScore(), 1));
        row.addCell(Toolbox.d2s(H.getScoreDelta(), 2));
        row.addCell(H.getRank() + "");
        row.addCell(H.getVerdict().toString());
        row.addCell(H.isDecoy() ? "Y" : " ");
        row.addCell(Toolbox.d2s(pt.getX(), 0));
        row.addCell(Toolbox.d2s(1E6 * pt.getY() / H.getMzObs(), 1));
        row.addCell(H.getCharge() + "");
        row.addCell(Toolbox.d2s(H.getMzObs(), 3));
        row.addCell(H.getPrevAA());
        Td C = new Td();
        C.addContentWithoutEscaping(peptide);
        row.addCell(C);
        row.addCell(H.getNextAA());
        row.addCell(H.getPeptideHitId().toString());
        row.addCell(Toolbox.d2s(H.getDetectibility(), 0));
        row.addCell(H.getDeNovoCandidate());
      }
    }
    report.addDiv(new Div(line));
    for (String lib : MsMsDictionary.sLibrary.keySet())
      report.addDiv(new Div(lib + ": " + MsMsDictionary.sLibrary.cells(lib).toString()));

    return report;
  }
  public XYPoint calcResidual(float devi, int charge)
  {
    float mass = devi * (float )charge;
    Range<Float> iso = mSettings.getPrecursorMzRange();
    XYPoint       pt = null;
    if (Toolbox.isSet(mPTMs) && !mSettings.getPrecursorMzRange().isEnclosed(devi))
    {
      Collection<PTM> slice = mPTMs.subset(mass + iso.getLower(), mass + iso.getUpper());
      if (Toolbox.isSet(slice))
      {
        pt = new SimplePoint(Double.MAX_VALUE, 0);
        for (PTM p : slice)
          if (Math.abs(p.getMonoDeltaMass() - mass) < pt.getX())
          { pt.setX(Math.abs(p.getMonoDeltaMass() - mass)); pt.setY(p.getMonoDeltaMass()); }

        pt.setX(Math.round(pt.getY())); pt.setY(devi - pt.getY() / (float )charge);
      }
    }
    if (pt == null)
    {
      pt = new SimplePoint();
      pt.setX(Math.round(mass)); pt.setY(devi - pt.getX() / (float )charge);
    }
    return pt;
  }
  class Analysis
  {
    float                       iDecoyOffset      = 0f;
    Map<Long, Float>            iPrecursorMatches;
    Multimap<Long, MsIon>       iFragmentMatches;
    SetMultimap<Long, Float> iMatches;
    Collection<SimpleCandidate> iCandidates;

    Analysis()             { super(); }
    Analysis(float offset) { super(); iDecoyOffset = offset; }

    synchronized public void run()
    {
      // set a clean slate
      iPrecursorMatches = new HashMap<Long, Float>();
      iFragmentMatches  = HashMultimap.create();
      iMatches          = TreeMultimap.create();
      iCandidates       = new TreeSet<SimpleCandidate>(SimpleCandidate.SORT_SCORE_KEY_DESEND);

      queryByPrecursors();
      // non-sense to place the min on the match per fragment?
      MultiTreeMap<Long, Float> matches = MsMsDictionary.cells(mLibName).queryByFragments(
        getSettings(), iFragmentMatches, iPrecursorMatches, mFrags.values(),
        iDecoyOffset, 0, Toolbox.isSet(mPrecursorMzs));

      yell("\n" + iFragmentMatches.size()  + " fragments found" + (iDecoyOffset > 0 ? " with decoy offset of " + iDecoyOffset : ""));
      yell(Toolbox.toString(mFrags.values(), ","));

      if (Toolbox.isSet(matches))
        for (Long id : matches.keySet())
          iMatches.putAll(id, matches.cells(id));

      if (Toolbox.isSet(getPtmpMap()) && Toolbox.isSet(iPrecursorMatches)) queryByMods();      // score with expanded ions!!
      survey(); // score the simple matches first
    }
    public void queryByPrecursors()
    {
      if (!Toolbox.isSet(getPrecursorMzs())) return;

      for (Integer charge : getPrecursorMzs().keySet())
      {
        Collection<Range<Float>> slices = enumeratePrecursor(getSettings(), getPtmMs(), getPrecursorMzs().cells(charge), charge == 0 ? 2 : charge);
        if (Toolbox.isSet(slices))
          for (Range<Float> slice : slices)
            iPrecursorMatches.putAll(MsMsDictionary.cells(mLibName).queryByPrecursor(slice.getLower() + iDecoyOffset, slice.getUpper() + iDecoyOffset, charge));
      }
      // for debugging
      //SpectralLib_Util.print(mDictionary, mPrecursors, "Precursors @" + Toolbox.toString(mPrecursorMzs, ","));
      yell("Precursors @" + Toolbox.toString(mPrecursorMzs.values(), ",") + (iDecoyOffset > 0 ? "+" + iDecoyOffset : "") + ": " + iPrecursorMatches.size() + " found.");
    }
    protected void queryByMods()
    {
      if (!Toolbox.isSet(mPrecursorMzs) || !Toolbox.isSet(mPTMs)) return;

      // go through the precrusor matches
      Range<Float>    isolation = mSettings.getPrecursorMzRange();
      //Collection<MsIon> shifted = new ArrayList<MsIon>();
      for (Long precursor : iPrecursorMatches.keySet())
      {
        Float delta = (float )mSettings.getEffectivePrecursorTol((double )iPrecursorMatches.cells(precursor));
        //if (iPrecursorMatches.cells(precursor) > 597.6 && iPrecursorMatches.cells(precursor) < 597.605)
        //{
        //  System.out.println();
        //}
        for (MsIon ion : getPrecursorMzs().values())
        {
          if (ion.getCharge() == null || ion.getCharge() == 0) continue;

          Float offset = (float )ion.getMz() - iPrecursorMatches.cells(precursor);
          if (!isolation.isEnclosed(offset))
          {
            // pull the plausible PTMs
            Collection<PTM> trials = getPtmpMap().subset(ion.getCharge() * (offset - delta), ion.getCharge() * (offset + delta));
            if (Toolbox.isSet(trials))
              for (PTM ptm : trials)
              {
                //if (ptm.getName().indexOf("Oxidation") >= 0)
                //{
                //  System.out.println();
                //}
                // we used to have 1 min match per query fragment, why, 20120503
                MultiTreeMap<Long, Float> M = MsMsDictionary.cells(mLibName).queryByFragments(
                  getSettings(), iFragmentMatches, iPrecursorMatches, mFrags.values(),
                  ptm.getMonoDeltaMass().floatValue() + iDecoyOffset, 0, Toolbox.isSet(mPrecursorMzs));
                for (Long id : M.keySet())
                  for (Float m : M.cells(id)) iMatches.put(id, m * -1f);
              }
          }
        }
      }
    }
    synchronized public void survey()
    {
      if (iMatches == null) return;

      Collection<Long> ids = new ArrayList<Long>(iMatches.keySet());
      for (Long id : ids)
      {
        SimpleCandidate match = new SimpleCandidate();
        float multiplier = 1f;

        match.setPrimaryKey(id);
        for (MsIon M : iFragmentMatches.cells(id))
          if (M.getMz() < 0) { multiplier = 2f; break; }

        int nmatch = iFragmentMatches.cells(id).size(),
          nfrag = (int )(mFrags.size() * multiplier),
          nsig = MsMsDictionary.cells(mLibName).getSignatureSize(id);

        match.setScore(      -1d * hypergeometric_oneSidedPValue(nmatch, nfrag, nsig, mBinSize));
        if (nsig > 1) nsig--; else if (nfrag > 1) nfrag--;
        match.setScoreDelta(-1d * hypergeometric_oneSidedPValue(nmatch, nfrag, nsig, mBinSize) - match.getScore());
        match.setRank(nsig);

        iCandidates.add(match);
      }
    }
    synchronized public Map<PeptideHit, MsMsAssignment> score(Map<PeptideHit, MsMsAssignment> candidates, Long id, boolean decoy)
    {
      if (iMatches == null) return candidates;

      MsMsAssignment assignment = MsMsDictionary.cells(mLibName).cells(id);
      PeptideHit          match = assignment.getAssignment().clone();

      match.setPeptideHitId(id); match.isDecoy(decoy);
      try
      {
        match.setRank(iMatches.cells(id).size());
      }
      catch (Exception e)
      {
        e.printStackTrace();
        throw new RuntimeException("Failed to set the match rank!", e);
      }
      match.setMzObs(assignment.getMsMs().getPrecursorMz());
      if (Toolbox.isSet(match.getNextAA()))
        match.setNextAA(match.getNextAA().split("/")[0]);

      if (Toolbox.isSet(mPrecursorMzs))
      {
        double delta = 0d, best = Double.MAX_VALUE;
        for (MsIon ion : mPrecursorMzs.values())
          if (Math.abs(ion.getMz() - match.getMzObs()) < best) { delta = ion.getMz() - match.getMzObs(); best = Math.abs(delta); }
        match.setMassDeviation(delta);
      }

      // figure out which of the fragments were matched
      Collection<MsIon> signatures = extractSignatureFragments(assignment, 8d, 50d, 7);
      String unmatched = null;
      for (MsIon F : signatures)
      {
        double    tol = mSettings.getEffectiveFragmentTol(F.getMz());
        boolean found = false;
        for (MsIon M : iFragmentMatches.cells(id))
        {
          if (Math.abs(Math.abs(M.getMz()) - F.getMz() + iDecoyOffset) <= tol) { found = true; break; }
        }
        if (!found) unmatched = Toolbox.extend(unmatched, Toolbox.d2s(F.getMz(), 1), ";");
      }
      float multiplier = 1f;
      for (MsIon M : iFragmentMatches.cells(id))
      {
        if (M.getMz() < 0) { multiplier = 2f; break; }
      }

      int nmatch = iFragmentMatches.cells(id).size(), nfrag = (int )(mFrags.size() * multiplier), nsig = signatures.size();
      match.setScore(      -1d * hypergeometric_oneSidedPValue(nmatch, nfrag, nsig, mBinSize));
      if (nsig > 1) nsig--; else if (nfrag > 1) nfrag--;
      match.setErrorPctVar(-1d * hypergeometric_oneSidedPValue(nmatch, nfrag, nsig, mBinSize) - match.getScore());
      match.setDetectibility((double )signatures.size());
      match.setDeNovoCandidate(unmatched);

      //if (match.getScore() < 4)
      //{
      //  System.out.println();
      //}
      candidates.put(match, assignment);
      return candidates;
    }
  }
  private static void yell(String s)
  {
    if (sVerbose) System.out.println(s); else sMessage.append(s + "\n");
  }
  public static List<Peak> extractSignatureFragments(PeakList msms, Double min_snr, double half_width, int tops)
  {
    if (msms == null || !Toolbox.isSet(msms.getRawData())) return null;

    double baseline = msms.getIntensities()XYPoint_Util.getMinY(msms.getRawData(), 0d, 10000d),
      min_merge_count = msms.getNumMergedSpectra() * 0.5d;
    List<Peak> signature = new ArrayList<Peak>(), orphans = new ArrayList<Peak>();
    for (ScorableMsPeak ion : (List< ScorableMsPeak> )msms.getData(new ScorableMsPeakFactory()))
    {
      // apply the deviation so we have the theoretical m/z
      if (ion.getMzDeviation() != null) ion.setMz(ion.getMz() - ion.getMzDeviation());
      if (MsIon_Util.isSignature(ion, 450d, msms.getPrecursorMz()))
      {
        double pk_counts = XYPoint_Util.countValid(msms.getRawData(), ion.getMz() - half_width, ion.getMz() + half_width),
          //base = XYPoint_Util.getBaseline(msms.getRawData(), ion.getMz() - 75d, ion.getMz() + 75d, (int )Math.min(7, pk_counts * 0.10), min_snr),
          base = XYPoint_Util.getBaseline(msms.getRawData(), ion.getMz() - half_width, ion.getMz() + half_width, tops, min_snr),
          //noise = pk_counts > 7 ? base : baseline,
          noise = Math.max(base, baseline * min_snr), // base is negative if less than tops peaks in the region
          snr = ion.getIntensity() / noise;
        if ((min_snr == null || ion.getIntensity() > noise) &&
          (ion.getMergeCountRaw() == null ||
            ion.getMergeCountRaw().intValue() == 0 ||
            ion.getMergeCountRaw().intValue() >= min_merge_count))
        {
          //ion.setIntensityBias((float )noise); ion.setHTML("neighbors: " + pk_counts); // for debugging only
          ion.setIntensityBias((float )Math.abs(base));
          signature.add(ion);
        }
        else if (pk_counts < 5) orphans.add(ion);
      }
    }
    if (signature.size() < 5) signature.addAll(orphans);
    return signature;
  }
  public static Multimap<String, MsMsAssignment> readNistPeptideLib(Multimap<String, MsMsAssignment> assignments,
    Date create_date,
    NistPepLibFile data_source)
    throws Exception
  {
    if (assignments == null) assignments = HashMultimap.create();

    MsMsAssignment assign = null;
    int            counts = 0;
    while (data_source != null && ((assign = data_source.newMsMsAssignment()) != null))
    {
      if (++counts % 2000   == 0) System.out.print(".");
      if (  counts % 250000 == 0) System.out.println(counts);

      assign.getAssignment().setProperty("CreationDate", create_date != null ? create_date : new Date());

      // only deposit the new one
      boolean found = false;
      if (Toolbox.isSet(assignments.cells(assign.getAssignment().getSequence())))
        for (MsMsAssignment as : assignments.cells(assign.getAssignment().getSequence()))
        {
          if (       assign.getAssignment().getCharge() == as.getAssignment().getCharge() &&
            Math.abs(assign.getMsMs().getPrecursorMz()   - as.getMsMs().getPrecursorMz()) < 0.5)
          {
            found = true; break;
          }
        }
      if (!found) assignments.put(assign.getAssignment().getSequence(), assign);
    }
    data_source.close();

    return assignments;
  }
  public static MultiTreeMap<String, SQUID_Sequence> mapNistPeptideLib2Genes(Collection<MsMsAssignment> assignments) throws Exception
  {
    Collection<String> accs = new HashSet<String>(), amgen_geneids = new HashSet<String>();
    for (MsMsAssignment as : assignments)
    {
      if      (as.getAssignment().getId().indexOf("sp") == 0) accs.add("swissprot_aa//" + as.getAssignment().getId());
      else if (as.getAssignment().getId().indexOf("tr") == 0) accs.add(   "trembl_aa//" + as.getAssignment().getId());
    }

    Retriever retriever = new Retriever();
    retriever.setSuppressDuplicates(true);

    List<SQUID_Sequence> sequences = retriever.batchRetrieval(accs);

    MultiTreeMap<String, SQUID_Sequence> acc_seq = new MultiTreeMap<String, SQUID_Sequence>();
    for (SQUID_Sequence s : sequences)
    {
      amgen_geneids.add(s.getAmgenGeneId());
      for (Accession acc : s.getAccessions())
        acc_seq.add(acc.getID(), s);
    }

    System.out.println(amgen_geneids.size() + " genes detected.");

    Map<String, GeneFeature> geneid_feature = fetchGeneDefinition(amgen_geneids);
    // populate the gene info for the sequences
    for (SQUID_Sequence seq : acc_seq.values())
      if (geneid_feature.cells(seq.getAmgenGeneId()) != null)
        seq.setProperty("GeneFeature", geneid_feature.cells(seq.getAmgenGeneId()));

    return acc_seq;
  }
  public static Map<String, GeneFeature> fetchGeneDefinition(Collection<String> geneids) throws SQLException
  {
    Connection conn = SQUIDConnectionMgr.getInstance().getConnectionFromPool(ConnectionMode.READONLY);
    List<List<String>> listOfLists = CollectionUtil.createSubLists(geneids, 250);

    Map<String, GeneFeature> geneid_feature = new HashMap<String, GeneFeature>();
    for (List<String> sublist : listOfLists)
    {
      SQLQuery query = new SQLQuery();
      query.addSelect("gene_id");
      query.addSelect("gene_symbol");
      query.addSelect("gene_name");
      query.addSelect("gene_description");
      query.addSelect("gene_type_name");
      query.addSelect("chromosome");
      query.addFromTable("SQUID.GENE_MV");
      query.addWhereClause("gene_id in (" + SQLUtil.generateSQLStringList(sublist) + ")");

      ResultSet rs = null;
      try
      {
        rs = query.execute(conn);
        while (rs.next())
        {
          GeneFeature gene = new GeneFeature(rs.getString(1), rs.getString(2), rs.getString(3)
          ).setDescription(rs.getString(4)).setGeneType(rs.getString(5)).setChromosome(rs.getString(6));
          geneid_feature.put(gene.getGeneId(), gene);
        }
      }
      finally { SQLUtil.closeResources(rs); }
    }
    return geneid_feature;
  }
  public static void doImportNistPeptideLib(Connection conn, Collection<PTM> ptms,
    Map<String, String> idmapping, Boolean save_by_src, Long expt_id, String datafile) throws Exception
  {
    NistPepLibFile source = new NistPepLibFile(datafile, ptms);
    Multimap<String, MsMsAssignment> assignments = HashMultimap.create();
    Multimap<String, String>            unmapped = HashMultimap.create();

    MsMsAssignment assign = null;
    int            counts = 0;
    while (source != null && ((assign = source.newMsMsAssignment()) != null))
    {
      if (idmapping != null)
      {
        String gene = null;
        for (String id : assign.getAssignment().getId().split("\\|"))
        {
          gene = idmapping.cells(id);
          if (gene != null) break;
        }
        if (gene == null && assign.getAssignment().getDescription().toUpperCase().indexOf("IG") >= 0)
        {
          gene = "IG_LIKE";
        }
        if (gene == null)
        {
          unmapped.put(assign.getAssignment().getId(), assign.getAssignment().getSequence());
        }
        assign.getAssignment().setProperty("Accession", assign.getAssignment().getId());
        assign.getAssignment().setId(gene);
      }

      if (assignments.size() < 5) assignments.put(assign.getMsMs().getRunTitle(), assign);

      if (++counts % 2000   == 0)
      {
        System.out.print(".");
        if (save_by_src != null && save_by_src)
        {
          for (String sample : assignments.keySet())
          {
            if (assignments.cells(sample).size() > 100)
            {
              Db_Util.batchInsert(conn, 123L, assignments.cells(sample));
              assignments.removeAll(sample);
            }
          }
        }
      }
      if (  counts % 250000 == 0) System.out.println(counts);
    }
    if (save_by_src != null && !save_by_src)
    {
      String run_title = new File(datafile).getName();
      Spectre_AnalysisRun run = Db_Util.newAnalysisRun(conn, expt_id, run_title.replaceAll("\\.msp", ""), "system", ProcessingStatus.CONFIGURING);
      Db_Util.batchInsert(conn, run.getAnalysisRunId(), assignments.values());
    }

    System.out.println(datafile);
    System.out.println("Unmapped Accessions:");
    for (String id : unmapped.keySet())
      System.out.println("|" + id + "|" + unmapped.cells(id).size() + "|");

    System.out.println("Break down by the instrument type:");
    for (MsInstrumentType inst : source.getInstrumentType().keySet())
      System.out.println(inst.toString() + ": " + source.getInstrumentType().cells(inst));

    System.out.println("\nBreak down by the Sample:");
    for (String sample : source.getSamples().keySet())
      System.out.println(sample + ": " + source.getSamples().cells(sample));

    source.close();
  }
  public static MapOfMultiMap<String, String, MsMsAssignment> doImportNistLib2011(Collection<PTM> ptms, Date creation_date, String datafile) throws Exception
  {
    Multimap<String, MsMsAssignment> assignments = readNistPeptideLib(null, creation_date, new NistPepLibFile(datafile, ptms));

    System.out.println("Assignment readSpectrumIdentifier: " + assignments.size());

    MultiTreeMap<String, SQUID_Sequence>                    acc_seq = mapNistPeptideLib2Genes(assignments.values());
    Multimap<String, MsMsAssignment>                     acc_assign = HashMultimap.create();
    MapOfMultiMap<String, String, MsMsAssignment> chrom_gene_assign = new MapOfMultiMap<String, String, MsMsAssignment>();

    for (MsMsAssignment as : assignments.values()) acc_assign.put(as.getAssignment().getId(), as);

    Set<SQUID_Sequence> seqs = new HashSet<SQUID_Sequence>();
    for (String accs : acc_assign.keySet())
    {
      seqs.clear();
      for (String acc : accs.split("\\|"))
        Toolbox.addNotNull(seqs, acc_seq.cells(acc));

      for (SQUID_Sequence seq : seqs)
      {
        GeneFeature gene = (GeneFeature )seq.getProperty("GeneFeature");
        for (MsMsAssignment A : acc_assign.cells(accs))
        {
          A.getAssignment().setProperty("ProtName", gene.getGeneSymbol()); // to become Spectre_PeptideAssignment.ProteinName later
          A.getAssignment().setDescription(gene.getGeneSymbol() + ": " + A.getAssignment().getDescription());
          chrom_gene_assign.add(gene.getChromosome(), gene.getGeneSymbol(), A);
        }
      }
    }
    return chrom_gene_assign;
  }
  public static void save(Connection conn, Long expt_id, MapOfMultiMap<String, String, MsMsAssignment> chrom_gene_assign) throws Exception
  {
    Spectre_Experiment parent = Spectre_Experiment.row(conn, expt_id);
    for (String ch : chrom_gene_assign.keySet())
    {
      Spectre_Experiment folder_ch = Db_Util.getExperiment(conn, parent, "ch " + ch);
      // go by the genes
      for (String gene : chrom_gene_assign.cells(ch).keySet())
      {
        Spectre_Experiment folder_gene = Db_Util.getExperiment(conn, folder_ch, gene.substring(0, 1) + "...");
        Spectre_AnalysisRun run = Db_Util.newAnalysisRun(conn, folder_gene.getRowId(), gene, "system", ProcessingStatus.CONFIGURING);
        Collection<MsMsAssignment> temp_as = new ArrayList<MsMsAssignment>();
        temp_as.add(Toolbox.front(chrom_gene_assign.cells(ch, gene)));
        //Db_Util.batchInsert(conn, run.getAnalysisRunId(), chrom_gene_assign.cells(ch, gene));
        Db_Util.batchInsert(conn, run.getAnalysisRunId(), temp_as);
      }
    }
  }
  public static Long saveBinaryLibs(DataOutput output, Collection<String> sequences, boolean compact, PeptideLibFile_Abstract... data_sources) throws Exception
  {
    if (!Toolbox.isSet(data_sources)) return 0L;

    MsMsAssignment assign = null;
    Long           counts = 0L;

    for (PeptideLibFile_Abstract lib : data_sources)
    {
      while (lib != null && ((assign = lib.newMsMsAssignment()) != null))
      {
        if (++counts % 5000   == 0) System.out.print(".");
        if (  counts % 500000 == 0) System.out.println(counts);

        //if (Toolbox.equals(assign.getAssignment().getSequence(), "YSVFFQSLAVVEQEMK"))
        //{
        //  System.out.println();
        //}
        if (assign == null || assign.getAssignment() == null) break;
        if (Toolbox.equals(assign.getStatus(), Subject.eStatus.skipped)) continue;

        String tag = assign.getAssignment().getBackboneCutPattern(false) + "z" + assign.getAssignment().getCharge();
        if (!sequences.contains(tag))
        {
          // empty the un-neccessary fields to save space
          assign.getMsMs().setRunTitle(null);
          assign.getMsMs().setScans((List<String> )null);

          // save the record
          if (compact) assign.compact_write(output); else assign.write(output);
          sequences.add(tag);
        }
        else
        {
          //System.out.println(tag);
        }
        assign.dispose();
      }
      lib.close();
    }
    return counts;
  }
  public static String Src2Bin(String binlib, PeptideLibFile_Abstract... msps) throws Exception
  {
    File lib = new File(binlib); if (lib != null && lib.exists()) lib.delete();
    BufferedRandomAccessFile  os = new BufferedRandomAccessFile(lib, "rw");
    Collection<String> sequences = new HashSet<String>();

    for (PeptideLibFile_Abstract msp : msps)
    {
      Long size = SpectralLib_Util.saveBinaryLibs(os, sequences, Toolbox.hasSubstr(binlib, ".compact"), msp);
      System.out.println("\n" + sequences.size() + "/" + size + ": " + msp.getName());
    }
    os.close();

    return lib.toString();
  }
  // prepare a BIN repository from experiments with peptide assignments
  public static Long saveBinaryLibs(DataOutput output, Collection<String> sequences, boolean compact, boolean consolidate,
    Connection conn, Collection<Long> runids) throws Exception
  {
    Collection<Long> msids = new HashSet<Long>();

    for (Long runid : runids)
      msids.addAll(Spectre_MsMsSpectrum.getRowIdsAnnotatedForAnalysisRun(conn, runid, 2, 5d, null));

    Multimap<String, MsMsAssignment> seq_assignment =
      Qualitative_Util.toMsMsAssignments(Spectre_MsMsSpectrum.getRowsAnnotated(conn, msids, MsMsSpectrumRetrievalFlags.ASSIGNMENT));

    for (String seq : seq_assignment.keySet())
    {
      // pick the most intense one by the BPI
      List<MsMsAssignment> pool = new ArrayList<MsMsAssignment>(seq_assignment.cells(seq));
      Collections.sort(pool, MsMsAssignment.BPI_DESCEND);

      if (compact) Toolbox.front(pool).compact_write(output); else Toolbox.front(pool).write(output);
    }
    return 0L;
  }
  public static void saveProteinIdsAsBinaryLib(Connection conn, String libfilename, Float mztol, ProteinMatrix_View matrix) throws Exception
  {
    System.out.println("Reading the protein IDs from the protein matrix:");
    System.out.println(matrix.getProteinIdReportFile().getName());
    System.out.println("and save the binary library to ");
    System.out.println(libfilename);

    //ProteinMatrix_View matrix = ProteinId_Util.newProteinMatrix_View("1331824947724.bin", "wyu");
    matrix.load(null);
    BufferedRandomAccessFile os = new BufferedRandomAccessFile(libfilename + ".bin", "rw");
    Collection<Long>       pids = new HashSet<Long>();

    int order = 0, assignments = 0, msms = 0;
    for (ProteinId pid : matrix.getCaller().getProteinIdGroup().keySet())
    {
      if (++order % 50 == 0) System.out.print(".");
      ProteinId sum = Qualitative_Util.toSummation(pid, matrix.getCaller().getAssignmentMap());
      if (sum.getAssignments().keySet().size() > 1)
      {
        for (String peptide : sum.getAssignments().keySet())
        {
          pids.clear();
          for (MsMsAssignment A : sum.getAssignments().cells(peptide)) pids.add(A.getId());

          Collection<Spectre_MsMsSpectrum> rows = Spectre_MsMsSpectrum.getRowsForPeptideAssignmentIds(conn, pids, true);
          rows = Spectre_MsMsSpectrum.populateSpectra(  conn, rows, true);
          rows = Spectre_MsMsSpectrum.attachAssignments(conn, rows, false, true, true);
          // remove the noise peaks using localDenoise. Merge count will be considered during signature extraction
          Spectre_MsMsSpectrum composit = SpectralClustering_Util.composite(rows, null, mztol, true);
          composit.getPeptideAssignment().matchToSpectrum(conn, composit.getRawSpectrum(), 0.5f, true);

          MsMsAssignment assignment = Toolbox.front(sum.getAssignments().cells(peptide));
          assignment.setMsMs(composit.convertToMsMsMatches());
          assignment.getMsMs().setScanType(MsScanType.CENTROID);
          assignment.write(os);
          assignments++;
          msms += sum.getAssignments().cells(peptide).size();
        }
      }
      //break; // for debugging
    }
    os.close();
    System.out.println("Binary Lib save (protein id/distinct assignments/msms:" + order + "/" + assignments + "/" + msms);
  }
  public static void saveProteinIdsAsBinaryLib(Connection conn, ProteinMatrix_View matrix, int row_limit, BinPepLibFile... libs) throws Exception
  {
    System.out.println("Reading the protein IDs from the protein matrix:");
    System.out.println(matrix.getProteinIdReportFile().getName());
    System.out.println("and save the binary library to ");
    System.out.println(Toolbox.toString(libs, ";\n"));

    matrix.load(null);
    saveProteinIdsAsBinaryLib(conn, matrix.getCaller(), null, null, row_limit, libs);
  }
  public static void saveProteinIdsAsBinaryLib(Connection conn, ProteinIDCaller caller,
    RandomAccessFile ds, Map<String, Long> peptide_pointer,
    int row_limit, BinPepLibFile... libs) throws Exception
  {
    saveProteinIdsAsBinaryLib(conn, caller.getProteinIdGroup().keySet(),
      caller.getAssignmentMap(), ds, peptide_pointer, row_limit, libs);

  }
  public static void saveProteinIdsAsBinaryLib(Connection conn, Collection<ProteinId> protein_ids,
    Map<String, MultiTreeMap<String, MsMsAssignment>> assignment_map,
    RandomAccessFile ds, Map<String, Long> peptide_pointer,
    int row_limit, BinPepLibFile... libs) throws Exception
  {
    int qualified = 0;
    for (ProteinId pid : protein_ids)
      if (pid.getAssignments().keySet().size() > 1) qualified++;

    System.out.println("Loading the assignments from " + qualified + " qualified protein IDs based on " +
      peptide_pointer.keySet().size() + " distinct peptides (pid/assignments/msms_merged/msms)");

    Collection<Long> pids = new HashSet<Long>();
    int order = 0, assignments = 0, msms = 0, msms_merged = 0;
    for (ProteinId pid : protein_ids)
    {
      if (++order % 50   == 0) System.out.print(".");
      if (  order % 1000 == 0) System.out.print(order + "/" + assignments + "/" + msms_merged + "/" + msms + "\n");
      //ProteinId sum = pid;
      if (pid.getAssignments().keySet().size() > 1)
      {
        ProteinId sum = ds == null ? Qualitative_Util.toSummation(pid, assignment_map) :
          Qualitative_Util.toSummation(pid, peptide_pointer, ds);
        for (String peptide : sum.getAssignments().keySet())
        {
          pids.clear();
          // trim the set by the error rate first
          for (MsMsAssignment A : sum.getAssignments().cells(peptide)) pids.add(A.getId());

          if (pids.size() > row_limit)
          {
            // no need to set the owner
            Map<Long, Spectre_MsMsSpectrum> pid_ms = Spectre_MsMsSpectrum.getRowMapWithTypesForPeptideAssignmentIds(conn, pids);
            for (MsMsAssignment A : sum.getAssignments().cells(peptide))
            {
              if (!pid_ms.containsKey(A.getId())) throw new RuntimeException("MS not found: " + A.getId());
              pid_ms.cells(A.getId()).setRetentionTime((float )A.getAssignment().getErrorPct());
            }
            List<Spectre_MsMsSpectrum> focused = new ArrayList<Spectre_MsMsSpectrum>(pid_ms.values());
            Collections.sort(focused, Spectre_MsMsSpectrum.RT_ASC);

            MapOfMultiMap<BinPepLibFile, Integer, Long> lib_z_pid = samplePeptideAssignments(focused, row_limit, libs);
            pids.clear(); pids.addAll(lib_z_pid.values());
          }

          Collection<Spectre_MsMsSpectrum> rows = Spectre_MsMsSpectrum.getRowsForPeptideAssignmentIds(conn, pids, false);
          rows = Spectre_MsMsSpectrum.populateSpectra(conn, rows, true);
          rows = Spectre_MsMsSpectrum.attachAttributes( conn, rows);
          rows = Spectre_MsMsSpectrum.attachAssignments(conn, rows, false, false, false);

          for (BinPepLibFile lib : libs)
            rows = savePeptideAssignments(conn, lib, rows, row_limit); // need no more than a dozen

          assignments++;
          msms += sum.getAssignments().cells(peptide).size();

          // dispose the temp objects
          Toolbox.dispose(rows);
        }
      }
      Toolbox.dispose(pid);
      //break; // for debugging
    }
    System.out.println("\n#entry\tthe Spectral Repositories");
    for (BinPepLibFile lib : libs)
    {
      System.out.println(lib.getNumberOfEntry() + ":\t\t" + lib.getName() + "(" + Toolbox.toString(lib.getInstrumentTypes(), ",") + "/" +
        Toolbox.toString(lib.getIonActivationTypes(), ",") + ")");
      lib.close();
    }

    System.out.println("\nBinary Lib save (protein id/qualifid/distinct assignments/msms:" + order + "/" + qualified + "/" + assignments + "/" + msms);
  }
  public static MapOfMultiMap<BinPepLibFile,Integer, Long>
  samplePeptideAssignments(Collection<Spectre_MsMsSpectrum> rows, int row_limit, BinPepLibFile... libs)
  {
    // remove the noise peaks using localDenoise. Merge count will be considered during signature extraction
    MapOfMultiMap<BinPepLibFile,Integer, Long> lib_charge_pid = new MapOfMultiMap<BinPepLibFile, Integer, Long>();
    for (Spectre_MsMsSpectrum ms : rows)
      for (BinPepLibFile lib : libs)
        if (lib_charge_pid.cells(lib, ms.getPrecursorCharge()) == null ||
          lib_charge_pid.cells(lib, ms.getPrecursorCharge()).size() < row_limit)
          if (lib.isCompatible(ms.getInstrumentType(), ms.getIonActivationType()))
            lib_charge_pid.add(lib, ms.getPrecursorCharge(), ms.getPeptideAssignmentId());

    return lib_charge_pid;
  }
  @Deprecated
  public static Collection<Spectre_MsMsSpectrum>
  savePeptideAssignments(Connection conn, BinPepLibFile lib,
    Collection<Spectre_MsMsSpectrum> rows, int row_limit) throws SQLException, IOException
  {
    if (!Toolbox.isSet(rows) || lib == null ||
      !Toolbox.isSet(lib.getInstrumentTypes()) ||
      !Toolbox.isSet(lib.getIonActivationTypes())) return rows;

    // remove the noise peaks using localDenoise. Merge count will be considered during signature extraction
    Multimap<Integer, Spectre_MsMsSpectrum> permitted = HashMultimap.create();
    Iterator<Spectre_MsMsSpectrum> itr = rows.iterator();
    while (itr.hasNext())
    {
      Spectre_MsMsSpectrum ms = itr.next();
      if (lib.isCompatible(ms.getInstrumentType(), ms.getIonActivationType()))
      {
        // always count on the charge state of the assignment
        ms.setPrecursorCharge(ms.getPeptideAssignment().getChargeState());
        if (row_limit <= 0 || permitted.cells(ms.getPrecursorCharge()) == null ||
          row_limit > permitted.cells(ms.getPrecursorCharge()).size())
          permitted.put(ms.getPrecursorCharge(), ms);
        // remove the spec from future considration
        itr.remove();
      }
    }
    // establish the average m/z tolerance from the instrument types
    Float mztol = 0f;
    for (MsInstrumentType inst : lib.getInstrumentTypes()) mztol += inst.getDefaultToleranceMs2();
    mztol /= lib.getInstrumentTypes().size();

    if (Toolbox.isSet(permitted))
    {
      //Spectre_MsMsSpectrum.populateSpectra(  conn, permitted.values(), true);
      //Spectre_MsMsSpectrum.attachAttributes( conn, permitted.values());
      //Spectre_MsMsSpectrum.attachAssignments(conn, permitted.values(), false, false, false);
      for (Integer z : permitted.keySet())
      {
        if (z == null || z == 0)
        {
          continue;
        }

        MsMsAssignment     assignment = Toolbox.front(permitted.cells(z)).toMsMsAssignment();
        Spectre_MsMsSpectrum composit = SpectralClustering_Util.composite(permitted.cells(z), null, mztol, false);
        composit.getPeptideAssignment().matchToSpectrum(conn, composit.getRawSpectrum(), mztol, true);

        assignment.setMsMs(composit.convertToMsMsMatches());
        assignment.getMsMs().setScanType(MsScanType.CENTROID);
        assignment.write(lib.getDataOutput());
        lib.increNumberOfEntry(1L);
      }
    }
    return rows;
  }
  public static void saveProteinIdsAsLibrary(Connection conn,float ms2tol, int sampling, String lib_root,
    Map<String, ProteinId> gene_pid,
    MapOfMultiMap<String, Integer, Long> peptide_z_asid,
    Table<String, String, MsMsAssignment> seq_peptidez_as)
    throws SQLException, IOException
  {
    Random rand = new Random(System.nanoTime());

    Table<MsInstrumentType, IonActivationType, BinPepLibFile> binlibs = newBinPepLibFilesMap(lib_root);
    MapOfMultiMap<MsInstrumentType, IonActivationType, Spectre_MsMsSpectrum> inst_act_ms =
      new MapOfMultiMap<MsInstrumentType, IonActivationType, Spectre_MsMsSpectrum>();

    System.out.println("\nSaving " + gene_pid.size() + " protein IDs to " + lib_root + " ...");

    int counts_pid = 0, as_totals = 0, nonqual_pids = 0;
    List<List<String>> genes = CollectionUtil.createSubLists(gene_pid.keySet(), 250);
    Set<String> abandoned = new HashSet<String>(), distinct = new HashSet<String>();
    for (List<String> gs : genes)
    {
      System.out.print(".");
      Collection<ProteinId> pids = Toolbox.getAll(gene_pid, gs);
      // batch the assignments
      Collection<Long> asids = new HashSet<Long>();
      for (ProteinId pid : pids)
      {
        pid.isInterest(true);
        if (pid.getRawAssignments().keySet().size() <= 1) { pid.isInterest(false); continue; }
        int promising = 0;
        if (pid.getRawAssignments().keySet().size() > 1)
          // have to look at the actual assignments, not those in the protein ID.
          for (String seq : pid.getRawAssignments().keySet())
            for (MsMsAssignment A : seq_peptidez_as.row(seq).values())
              if (A.getConcensusSize() > 2 && A.getAssignment().getErrorPct() < 5) promising++;

        if (promising == 0) { pid.isInterest(false); continue; }

        if (sampling > 0)
          for (String seq : pid.getRawAssignments().keySet())
            for (MsMsAssignment A : seq_peptidez_as.row(seq).values())
              for (Integer z : peptide_z_asid.cells(A.getAssignment().getBackboneCutPattern(false)).keySet())
                Toolbox.subsample(asids, peptide_z_asid.cells(A.getAssignment().getBackboneCutPattern(false), z), sampling, rand);
      }
      Map<Long, Spectre_MsMsSpectrum> id_spec = Spectre_MsMsSpectrum.getRowMapForPeptideAssignmentIds(conn, asids, false);

      if (sampling > 0)
        Spectre_MsMsSpectrum.populateSpectra(  conn, id_spec.values(), false);

      Spectre_MsMsSpectrum.populatePrecursorChargeFromAssignments(conn, id_spec.values());
      Spectre_MsMsSpectrum.attachAssignments(conn, id_spec.values(), false, true, true);

      for (ProteinId pid : pids)
        if (pid.isInterest())
        {
          counts_pid++;
          for (String seq : pid.getRawAssignments().keySet())
            if (!distinct.contains(seq))
              for (MsMsAssignment A : seq_peptidez_as.row(seq).values())
              {
                as_totals++; distinct.add(seq);
                // cells the protein info
                A.getAssignment().setDescription(pid.getProtein().getDescription());
                A.getAssignment().setId(         pid.getProtein().getId());
                MultiTreeMap<Integer, Long> z_id = peptide_z_asid.cells(A.getAssignment().getBackboneCutPattern(false));
                if (!Toolbox.isSet(z_id))
                {
                  System.out.println("Empty asid: " + A.toString());
                  continue;
                }
                inst_act_ms.clear();
                for (Long asid : z_id.values())
                {
                  Spectre_MsMsSpectrum ms = id_spec.cells(asid);
                  if (ms != null)
                    inst_act_ms.add(ms.getInstrumentType(), ms.getIonActivationType(), ms);
                }
                savePeptideAssignments(conn, inst_act_ms, binlibs, ms2tol, A);
              }
        }
        else
        {
          nonqual_pids++;
          abandoned.addAll(pid.getRawAssignments().keySet());
        }
      // displose the objects
      Toolbox.dispose(id_spec);
    }
    Toolbox.dispose(inst_act_ms);

    System.out.println("\nProtein IDs (qualified/nonqual/totals), Peptide Assignments/abandoned distinct peptides: " + counts_pid + "/" + nonqual_pids + "/" + gene_pid.size() + ", " + as_totals + "/" + abandoned.size());

    for (BinPepLibFile lib : binlibs.values())
    {
      if (lib.getDataOutput() != null) System.out.println(lib.getAbsolutePath() + " with " + lib.getNumberOfEntry() + " of entries created.");
      System.out.println("    # distinct (peptides/peptide_z/peptide seq): " + lib.mPeptides.size() + "/" + lib.mPeptideZs.size() + "/" + lib.mPeptideSeqs.size());
      lib.close();
    }
  }
  public static void savePeptideAssignments(Connection conn,
    MapOfMultiMap<MsInstrumentType, IonActivationType, Spectre_MsMsSpectrum> inst_act_ms,
    Table<MsInstrumentType, IonActivationType, BinPepLibFile> inst_act_libs,
    float ms2tol, MsMsAssignment A) throws IOException, SQLException
  {
    // save the assignments by the types
    for (MsInstrumentType inst : inst_act_ms.keySet())
      for (IonActivationType act : inst_act_ms.labelSet())
        if (inst_act_ms.cells(inst, act) != null)
        {
          // take the charge state from the assignment
          Multimap<Integer, Spectre_MsMsSpectrum> z_ms = Qualitative_Util.toChargeMap(inst_act_ms.cells(inst, act));
          for (Integer z : z_ms.keySet())
          {
            if (z == 0)
            {
              System.out.println("zero charge, " + A.toString());
              continue;
            }
            // set the charge state
            try
            {
              A.getAssignment().setCharge(z);
              Spectre_MsMsSpectrum composit = SpectralClustering_Util.composite(z_ms.cells(z), null, ms2tol, false);
              composit.setPeptideAssignment(Qualitative_Factory.newSpectrePeptideAssignment(A));
              composit.getPeptideAssignment().matchToSpectrum(conn, composit.getRawSpectrum(), ms2tol, true);

              if (z_ms.cells(z).size() > 2)
                MsIon_Util.purgeByMergeCount(composit.getRawSpectrum().getModifiableData(), 0.5d);

              A.setMsMs(composit.convertToMsMsMatches());
              A.getMsMs().setScanType(MsScanType.CENTROID);
              if (A.getAssignment().getProperties() != null)
                A.getAssignment().removeProperty("asid");
              if (inst_act_libs.cells(inst, act) != null)
              {
                A.write(inst_act_libs.cells(inst, act).getDataOutput());
                if (A.getMsMs() == null || !Toolbox.isSet(A.getMsMs().getData()))
                {
                  System.out.println("fragment ions not initialized?");
                }
                if (inst_act_libs.cells(inst, act).mPeptideZs.contains(A.getAssignment().getUniqueTag()))
                {
                  System.out.print("W?");
                }
                inst_act_libs.cells(inst, act).increNumberOfEntry(1L);
                inst_act_libs.cells(inst, act).mPeptides.add( A.getAssignment().getBackboneCutPattern(false));
                inst_act_libs.cells(inst, act).mPeptideZs.add(A.getAssignment().getUniqueTag());
                inst_act_libs.cells(inst, act).mPeptideSeqs.add( A.getAssignment().getSequence());
              }
              else throw new RuntimeException("Unable to write an assignment with " + inst + " and " + act + " to the lib output!");
            }
            catch (Exception ue) { ue.printStackTrace(); }
          }
        }
  }
  public static void saveSummaryProteinIds(Connection conn, String outroot, Long... expt_ids) throws Exception
  {
    ProteinIDCaller caller  = ProteinId_Util.newProteinIDCallerByExperiments(conn, null, expt_ids);
    caller.getSettings().toGeneratePeptideMatrix(false);

    System.out.println("Calling the protein matrix for expt#" + Toolbox.toString(expt_ids, ","));
    saveProteinIds(ProteinId_Util.doProteinID(conn, caller), outroot);
  }
  public static void saveProteinIds(ProteinIDCaller caller, String outroot) throws Exception
  {
    BufferedRandomAccessFile indices = new BufferedRandomAccessFile(outroot + "/single.pid",    "rw"),
      qualified = new BufferedRandomAccessFile(outroot + "/qualified.pid", "rw");

    int single_count = 0;
    //Map<String, Investigation> expts = new HashMap<String, Investigation>();
    //for (Investigation s : caller.getInvestigations()) expts.put(s.getName(), s);

    for (ProteinId pid : caller.getProteinIdGroup().keySet())
    {
      ProteinId sum = Qualitative_Util.toSummation(pid, caller.getAssignmentMap());
      if (sum.getAssignments().keySet().size() == 1)
      {
        single_count++;
        sum.write(indices);
      }
      else
      {
        sum.write(qualified);
      }
    }

    indices.close(); qualified.close();
    System.out.println("# of Protein ID (single/total): " + single_count + "/" + caller.getProteinIdGroup().keySet().size());
    System.out.println(outroot + "/single.pid");
    System.out.println(outroot + "/qualified.pid");
  }
  public static BinPepLibFile[] newBinPepLibFilesByInstIonActivation(String root)
  {
    BinPepLibFile it_cid = new BinPepLibFile(root + "_it_cid.bin").setInstrumentTypes(
      MsInstrumentType.FTLTQ, MsInstrumentType.LCQ, MsInstrumentType.LTQ, MsInstrumentType.ORBITRAP).setIonActivationTypes(IonActivationType.CID),
      qtof_cid = new BinPepLibFile(root + "_qtof_cid.bin").setInstrumentTypes(
        MsInstrumentType.QTOF, MsInstrumentType.ABI5600).setIonActivationTypes(IonActivationType.CID),
      it_etd = new BinPepLibFile(root + "_it_etd.bin").setInstrumentTypes(
        MsInstrumentType.FTLTQ, MsInstrumentType.LCQ, MsInstrumentType.LTQ, MsInstrumentType.ORBITRAP).setIonActivationTypes(IonActivationType.ETD, IonActivationType.ECD),
      it_hcd = new BinPepLibFile(root + "_it_hcd.bin").setInstrumentTypes(
        MsInstrumentType.FTLTQ, MsInstrumentType.LCQ, MsInstrumentType.LTQ, MsInstrumentType.ORBITRAP).setIonActivationTypes(IonActivationType.HCD);

    return new BinPepLibFile[] {it_cid, qtof_cid, it_etd, it_hcd};
  }
  private static Table<MsInstrumentType, IonActivationType, BinPepLibFile> put(Table<MsInstrumentType, IonActivationType, BinPepLibFile> tbl, BinPepLibFile... libs)
  {
    for (BinPepLibFile lib : libs)
      for (MsInstrumentType inst : lib.getInstrumentTypes())
        for (IonActivationType act : lib.getIonActivationTypes()) tbl.put(inst, act, lib);

    return tbl;
  }

  public static Table<MsInstrumentType, IonActivationType, BinPepLibFile> newBinPepLibFilesMap(String root)
  {
    Table<MsInstrumentType, IonActivationType, BinPepLibFile> tbl = HashBasedTable.create();

    return put(tbl, new BinPepLibFile(root + "_it_cid.bin").setInstrumentTypes(
        MsInstrumentType.FTLTQ, MsInstrumentType.LCQ, MsInstrumentType.LTQ, MsInstrumentType.ORBITRAP).setIonActivationTypes(IonActivationType.CID),
      new BinPepLibFile(root + "_qtof_cid.bin").setInstrumentTypes(
        MsInstrumentType.QTOF, MsInstrumentType.ABI5600).setIonActivationTypes(IonActivationType.CID),
      new BinPepLibFile(root + "_it_etd.bin").setInstrumentTypes(
        MsInstrumentType.FTLTQ, MsInstrumentType.LCQ, MsInstrumentType.LTQ, MsInstrumentType.ORBITRAP).setIonActivationTypes(IonActivationType.ETD, IonActivationType.ECD),
      new BinPepLibFile(root + "_it_hcd.bin").setInstrumentTypes(
        MsInstrumentType.FTLTQ, MsInstrumentType.LCQ, MsInstrumentType.LTQ, MsInstrumentType.ORBITRAP).setIonActivationTypes(IonActivationType.HCD));
  }
  // "/bioinfo/scratch/wyu/staging/sequences/cp_reviewed_9606.fasta"
  public static List<ProteinId> Bins2Bin(BinPepLibFile outfile, String dbfile, BinPepLibFile... bins) throws Exception
  {
    Set<String> peptides = new HashSet<String>();
    for (BinPepLibFile lib : bins)
    {
      System.out.println("Reading the assignments from " + lib.getName() + "...");
      int counts = 0;
      try
      {
        MsMsAssignment assign = lib.newMsMsAssignment();
        while (assign != null && assign.getMsMs() != null)
        {
          if (++counts % 50000 == 0) System.out.print(".");

          String peptide = assign.getAssignment().getUniqueTag();
          if (!peptides.contains(peptide))
          {
            assign.write(outfile.getDataOutput());
            peptides.add(peptide);
          }
          assign.dispose();
          assign = lib.newMsMsAssignment();
        }
      }
      catch (Exception e)
      {
        e.printStackTrace();
        throw e;
      }
      finally { lib.close(); }
      System.out.println(counts);
    }
    outfile.close();

    List<ProteinId> pids = null;
    if (Toolbox.isSet(dbfile))
    {
      System.out.println("Prep the protein IDs according to " + dbfile + "...");
      pids = ProteinId_Util.doProteinIds(Qualitative_Util.newProteinIDs(
        Sequence_Util.newProteins(dbfile)), outfile);
    }
    return pids;
  }
  public static void prepareMimslLib(String libname, String dbfile, double min_snr, double half_width, int ntops, BinPepLibFile... bins) throws Exception
  {
    BinPepLibFile combined = new BinPepLibFile(libname + ".bin");
    SpectralLib_Util.Bins2Bin(combined, dbfile, bins);

    MsMsDictionary dic = new MsMsDictionary().index(combined, min_snr, half_width, ntops);
    BufferedRandomAccessFile indices = new BufferedRandomAccessFile(libname + ".indices", "rw");

    dic.write(indices);
    indices.close();
  }
  public int getSignatureSize(Long id) { return mIdSigsize != null ? mIdSigsize.cells(id) : null; }
  public static MsMsDictionary cells(String name)
  {
    // settle for no PTM for now because we can;t cells to the Spectre_PTM outside of SPECTRE
    return prepare(name, true, null);
  }
  public static MsMsDictionary prepare(String name, boolean wait, Collection<PTM> ptms)
  {
    if (!Toolbox.isSet(name)) throw new RuntimeException("NULL library");

    MsMsDictionary dic = sLibrary.cells(name);
    if (dic == null)
    {
      // do we have it started already
      Future<MsMsDictionary> future = sFutures.cells(name);
      if (future == null)
        sFutures.put(name, Executors.newFixedThreadPool(1).submit(
          //new MsMsDictionary(name, Spectre_PTM.getAllasPTM().values())));
          new MsMsDictionary(name, ptms)));

      if (wait)
      {
        try
        {
          sLibrary.put(name, sFutures.cells(name).cells()); // let's wait for it
          sFutures.remove(name);
        }
        catch (Exception e) { throw new RuntimeException(e); }
      }
    }

    return sLibrary.cells(name);
  }
  public static boolean isReady(String name)
  {
    return (sFutures.cells(name) == null || sFutures.cells(name).isDone());
  }
  public long size() { return mPrecursors.size(); }
  public MsMsDictionary setName(String s) { mName = s; return this; }
  public String getName() { return mName; }
  public MsMsDictionary setSpectralLib(String s) { mSpectralLib = s; return this; }

  synchronized public MsMsAssignment cells(Long offset)
  {
    String binfile = libsdir + mSpectralLib;
    try
    {
      //System.out.println("assignment @" + offset);
      if (mBinLibrary == null)
      {
        File bin = new File(binfile);
        if (!bin.exists()) throw new RuntimeException("Spectra repository not found: " + binfile);

        mBinLibrary = new BufferedRandomAccessFile(binfile, "r");
      }
      mBinLibrary.seek(offset);

      return new MsMsAssignment(mBinLibrary);
    }
    catch (IOException e) { throw new RuntimeException("Not able to access the spectra repository: " + binfile, e); }
  }
  public void close()
  {
    try
    {
      if (mBinLibrary != null) { mBinLibrary.close(); mBinLibrary = null; }
    }
    catch (IOException e) { return; }
  }

  public MsMsDictionary add(Long id, MsIon precursor, Collection<MsIon> fragments, PeptideHit hit)
  {
    if (Toolbox.isSet(fragments))
    {
      if (mFragments  == null) mFragments  = new IonDictionary(2);
      if (mPrecursors == null) mPrecursors = new IonDictionary(3);
      if (mIdSigsize  == null) mIdSigsize  = new HashMap<Long, Integer>();
      //if (mChargeSequence == null) mChargeSequence = TreeMultimap.create();

      mPrecursors.add(precursor.getCharge(), (float )precursor.getMz(), id);
      mFragments.addAll(fragments, id);
      mIdSigsize.put(id, fragments.size());
    }

    return this;
  }

  //** Fetch the pool of peptide keys by the precursor m/z values. THe tolerance can be defined by ppm or isolation alone or both.
  public Map<Long, Float> queryByPrecursor(Float low, Float high)
  {
    return queryByPrecursor(low, high, 0);
  }
  public Map<Long, Float> queryByPrecursor(Float low, Float high, Integer z)
  {
    if (mPrecursors == null)     throw new RuntimeException("NULL  Precursor Index!");
    if (mPrecursors.size() == 0) throw new RuntimeException("Empty Precursor Index!");

    if (low == null || high == null) return null;
    return mPrecursors.query(low, high, z, null);
  }
  public MultiTreeMap<Long, Float> queryByFragments(MsSettings settings,
    Multimap<Long, MsIon> fragment_matches,
    Map<Long, Float> precursors, Collection<MsIon> frags,
    Float ptm, int min_count, boolean hasPrecursor)
  {
    if (mFragments == null)     throw new RuntimeException("NULL  Fragment Index!");
    if (mFragments.size() == 0) throw new RuntimeException("Empty Fragment Index!");

    if (!Toolbox.isSet(frags)) throw new RuntimeException("Query Fragments not initialized!");

    //    boolean              hasPrecursor = Toolbox.isSet(precursors);
    MultiTreeMap<Long, Float> matches = new MultiTreeMap<Long, Float>();
    for (MsIon ion : frags)
    {
      float                 delta = (float )settings.getEffectiveFragmentTol(ion.getMz());
      Map<Long, Float> candidates = mFragments.query((float )(ion.getMz() - delta), (float )(ion.getMz() + delta), ion.getCharge(), ptm);

      if (candidates.size() > min_count)
        for (Long id : candidates.keySet())
        {
          if (!hasPrecursor || precursors.keySet().contains(id)) matches.add(id, (float )ion.getMz());
          fragment_matches.put(id, ion);
        }
    }
    return matches;
  }
  public void print(Collection<Long> ids, String title)
  {
    System.out.println("\n" + title);
    for (Long id : ids)
    {
      MsMsAssignment assignment = cells(id);
      if (assignment != null)
        System.out.println("m/z" + Toolbox.d2s(assignment.getMsMs().getPrecursorMz(), 2) + ",+" +
          assignment.getMsMs().getPrecursorCharge() + ", " +
          assignment.getAssignment().getPeptideCutPattern());
    }
  }

  public String toString()
  {
    StringBuffer buf = new StringBuffer("MsMs Dictionary built from the following Spectral Libraries:\n" + mSpectralLib + "\n");

    buf.append("\n      Total Ms/Ms Spectra: " + mPrecursors.getIonID().size());
    buf.append("\nTotal Signature Fragments: " +  mFragments.getIonID().size());

    return buf.toString();
  }

  public MsMsDictionary init(String libname, Collection<PTM> ptms) throws Exception
  {
    long         msec0 = System.currentTimeMillis();
    long      freemem0 = Runtime.getRuntime().freeMemory();
    MsMsDictionary dic = new MsMsDictionary();

    String lib = sLibraryName.cells(libname);
    if      (lib == null) throw new RuntimeException("Unknown library: " + libname);

    File libfile = new File(libsdir + lib + ".indices");

    if (!libfile.exists()) throw new RuntimeException("The library index file was not found: " + libfile.getAbsolutePath());

    BufferedRandomAccessFile indices = new BufferedRandomAccessFile(libfile, "r"); // rw?
    dic.readSpectrumIdentifier(indices);
    if (dic != null) dic.setName(libname);

    indices.close();

    System.gc();
    System.out.println("     Elipse Time (sec) " + Toolbox.d2s((System.currentTimeMillis() - msec0) * 1E-3, 1));
    System.out.println("Free Memory Delta (MB) " + Toolbox.d2s((Runtime.getRuntime().freeMemory() - freemem0) * 1E-6, 1));

    System.out.println(dic.toString());
    sMessage.append(dic.toString());

    return dic;
  }

  public MsMsDictionary call() throws Exception
  {
    return init(mName, mPTMs);
  }
  public static MsMsDictionary index(String libname, double min_snr, double half_width, int tops, String... sequences) throws Exception
  {
    MsMsDictionary dic = new MsMsDictionary();
    dic.index(new BinPepLibFile(libsdir + sLibraryName.cells(libname) + ".bin"), min_snr, half_width, tops, sequences);

    return dic;
  }

  public static boolean isY(ScorableMsPeak ion)
  {
    return ion != null && (ion.is(MsDataFlag.FLAG_Y, MsDataFlag.FLAG_Y2, MsDataFlag.FLAG_Y_MINUS_NEUTRAL));
  }
  public static boolean isB(ScorableMsPeak ion)
  {
    return ion != null && (ion.is(MsDataFlag.FLAG_B, MsDataFlag.FLAG_B2, MsDataFlag.FLAG_B_MINUS_NEUTRAL));
  }
  public static boolean isCterm(ScorableMsPeak ion)
  {
    return isY(ion) || ion.is(MsDataFlag.FLAG_Y_MINUS_NEUTRAL, MsDataFlag.FLAG_Y_MINUS_NH3);
  }
  public static boolean isNterm(ScorableMsPeak ion)
  {
    return isB(ion) || (ion.is(MsDataFlag.FLAG_B_MINUS_H2O, MsDataFlag.FLAG_B_MINUS_NEUTRAL, MsDataFlag.FLAG_A));
  }
  public static boolean isYB(ScorableMsPeak ion)
  {
    //return ion != null && (ion.is(MsDataFlag.FLAG_Y, MsDataFlag.FLAG_Y2, MsDataFlag.FLAG_B, MsDataFlag.FLAG_B2));
    return isY(ion) || isB(ion);
  }
  public static boolean isAux(ScorableMsPeak ion)
  {
    return ion != null && ion.is(MsDataFlag.FLAG_IGNORE, MsDataFlag.FLAG_NOISE, MsDataFlag.FLAG_NOT_LABELED,
      MsDataFlag.FLAG_PRECURSOR, MsDataFlag.FLAG_PRECURSOR_NEUTRAL);
  }
  public static boolean isNeglectable(ScorableMsPeak ion)
  {
    return ion != null && ion.is(MsDataFlag.FLAG_IGNORE, MsDataFlag.FLAG_NOISE, MsDataFlag.FLAG_NOT_LABELED);
  }
  public static boolean isAnnotated(ScorableMsPeak ion)
  {
    return ion != null && Toolbox.isSet(ion.getLabel()) && !ion.getLabel().equals("?");
  }
  public static boolean isSignature(ScorableMsPeak ion, double min_mz, double precursor_mz)
  {
    double mass = ion.getCharge() != null && ion.getCharge() != 0 ? Mass.mzToMolWeight(ion.getMz(), ion.getCharge()) : ion.getMz();
    // some of the published transition use 2+ ions below precursor m/z
    return mass > precursor_mz && ion.getMz() > min_mz && !isAux(ion) && isYB(ion);
  }
  */
}
