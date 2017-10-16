package org.ms2ms.data.ms;

import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.Disposable;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.math.Histogram;
import org.ms2ms.math.QVals;
import org.ms2ms.math.Stats;
import org.ms2ms.mzjava.AnnotatedPeak;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

/** Keeper of the peptide matches to a MS/MS spectrum. It's not a single PSM since we don;t assume a single precursor
 *
 * Created by yuw on 8/7/16.
 */
public class Ms2Hits_ implements Disposable
{
  public static final String REJECT_PEAKSPLITTING = "rejected due to peak splitting";
  public static final String REJECT_SKEWED_ISO = "rejected due to skewed isotope envelop";
  public static final String REJECT_SPARSEPEAK    = "rejected due to insufficient peaks";
  public static final String REJECT_NONE          = "not rejected";
  public static final String CTR_GAP     = "GapScoreCenter";
  public static final String SIG_GAP     = "GapScoreSigma";

  private Map<String, ScoreModel> mScoreModels;
  private String mRejection;
  private Ms2Hit_ mTopRanked=null, mBestDecoy=null;
  private Map<String, Object> mPeakStats;
  private Map<String, Double> mBasis = new HashMap<>();
  private Map<String, Ms2Hit_> mExceptionals = new HashMap<>();

  private Histogram mDecoyY=null, mDecoyB=null;

  private MsnSpectrum mSpectrum;
  // 'ladders' build from the y and b ions, respectively (Score, DistinctSeqID, FpmEntry).
  private Table<Double, Integer, FpmEntry> mCtSegments=null, mNtSegments;
  private Multimap<Double, Ms2Hit_>   mExactMatches, mOpenMatches, mCandidates, mFinalists;

  public Ms2Hits_() { super(); }
  public Ms2Hits_(Map<String, Object> stats)
  {
    super(); mPeakStats=stats;
  }

  public Ms2Hits_(MsnSpectrum ms, Table<Double, Integer, FpmEntry> n, Table<Double, Integer, FpmEntry> c)
  {
    super();
    mSpectrum=ms; mCtSegments=c; mNtSegments=n; initFpmEntries();
  }

//  public boolean hasExceptional(String s) { return mExceptionals!=null && mExceptionals.get(s)!=null; }
  public Ms2Hit_ getExceptional(String s) { return mExceptionals.get(s); }
  public Map<String, Double> getScores() { return mBasis; }
  public boolean hasBasis(String... s)
  {
    if (!Tools.isSet(mBasis)) return false;
    if (Tools.isSet(s))
      for (String t : s) if (!mBasis.containsKey(t)) return false;

    return true;
  }
  public ScoreModel getScoreModel(String s) { return mScoreModels!=null?mScoreModels.get(s):null; }
//  public Map<String, Double> getBasis() { return mBasis; }
  public Double getBasis(String s) { return mBasis.get(s); }
  public double getGapScoreZ(double s)
  {
    return hasBasis(CTR_GAP, SIG_GAP) ? ((s-getBasis(CTR_GAP))/(getBasis(SIG_GAP)!=null?getBasis(SIG_GAP):1d)):s;
  }
//  public double getMatchProbZ(double s)
//  {
//    return hasBasis(CTR_MATCH, SIG_MATCH) ? ((s-getBasis(CTR_MATCH))/(getBasis(SIG_MATCH)!=null?getBasis(SIG_MATCH):1d)):s;
//  }
  public Map<String, Object> getPeakStats() { return mPeakStats; }

  public MsnSpectrum getSpectrum() { return mSpectrum; }

  public Table<Double, Integer, FpmEntry> getCtSegments() { return mCtSegments; }
  public Table<Double, Integer, FpmEntry> getNtSegments() { return mNtSegments; }

  public Multimap<Double, Ms2Hit_> getExactMatches() { return mExactMatches; }
  public Multimap<Double, Ms2Hit_> getOpenMatches()  { return mOpenMatches; }
  public Multimap<Double, Ms2Hit_> getCandidates()   { return mCandidates; }
  public Multimap<Double, Ms2Hit_> getFinalists()   { return mFinalists; }

  public Ms2Hit_ getBestDecoy() { return mBestDecoy; }

  public Collection<Ms2Hit_> getMatches()
  {
    Collection<Ms2Hit_> mm = new ArrayList<>(size());
    if (Tools.isSet(getExactMatches())) for (Ms2Hit_ H : getExactMatches().values()) mm.add(H);
    if (Tools.isSet(getOpenMatches()))  for (Ms2Hit_ H : getOpenMatches().values())  mm.add(H);

    return mm;
  }
  public Ms2Hits_ setSpectrum(MsnSpectrum s) { mSpectrum=s; return this; }

  public Ms2Hits_ setFinalists(Multimap<Double, Ms2Hit_> s)   { mFinalists=s; return this; }

  //  public Ms2Hits_ putExceptional(String s, Ms2Hit_ H) { mExceptionals.put(s, H); return this; }
  public Ms2Hits_ setScore(String k, Double s) { if (k!=null && s!=null) mBasis.put(k, s); return this; }
  public Ms2Hits_ setPeakCounts(Map<String, Object> s) { mPeakStats=s; return this; }
  public Ms2Hits_ reject(String reason) { mRejection=reason; return this; }
  public boolean isRejected() { return !(mRejection==null || mRejection.equals(REJECT_NONE)); }
  public String getRejection() { return mRejection; }

  public Ms2Hits_ setOpenMatches( Multimap<Double, Ms2Hit_> s) { mOpenMatches =null; mOpenMatches =s; return this; }
  public Ms2Hits_ setExactMatches(Multimap<Double, Ms2Hit_> s) { mExactMatches=null; mExactMatches=s; return this; }

  public Ms2Hits_ addMatch(Ms2Hit_ s, Tolerance tol)
  {
    return addMatch(s, tol.withinTolerance(s.getCalcMH(), s.getCalcMH()+s.getDelta()));
  }
  public Ms2Hits_ addMatch(Ms2Hit_ s, boolean isExact)
  {
    if (mExactMatches==null) mExactMatches = TreeMultimap.create(Ordering.natural().reverse(), Ordering.natural());
    if (mOpenMatches ==null) mOpenMatches  = TreeMultimap.create(Ordering.natural().reverse(), Ordering.natural());

    if (isExact) mExactMatches.put(s.getGapScore(), s.isExact(true));
    else          mOpenMatches.put(s.getGapScore(), s.isExact(false));

    return this;
  }
  public Ms2Hits_ addFinalist(Ms2Hit_ s)
  {
    if (mFinalists==null) mFinalists = TreeMultimap.create(Ordering.natural().reverse(), Ordering.natural());
    // sort the candidates by the e-val
    mFinalists.put(s.getScore(), s);

    return this;
  }

  public Multimap<Double, Ms2Hit_> consolidate(Multimap<Double, Ms2Hit_> hits, HashMap<Integer, Ms2Hit_> distincts, MultiTreeTable<Integer, String, Ms2Hit_> seqs,
                                         float[] AAs, OffsetPpmTolerance tol, Range<Integer> isoErr, float deci, TreeMultimap<Float, String> blocks)
  {
    // start from the Exact matches, higher Gap score first
    if (!Tools.isSet(hits)) return hits;

    Iterator<Map.Entry<Double, Ms2Hit_>> itr = hits.entries().iterator();
    while (itr.hasNext())
    {
      Map.Entry<Double, Ms2Hit_> E = itr.next();
      // potential matches
      List<Integer> hashes = E.getValue().hashcodeByIntervals(AAs, tol, isoErr, deci, blocks);
      Ms2Hit_ H = E.getValue();
      // check against the case where the same mass ladder is used on multiple peptides
      Integer mzHash = H.hashcodeByYBmz();
      hashes.add(mzHash);

      // need to find out the prior hit(s)
      Boolean redundant=null;
      for (Integer hash : hashes)
        if (distincts.containsKey(hash))
        {
          if (H.getSequence().equals(distincts.get(hash).getSequence()) &&
              H.getScore()          >distincts.get(hash).getScore()+5d ||
             (H.getProteinKey()>0 && distincts.get(hash).getProteinKey()<0)) // in case of small peptide
          {
            // more complicate if we need to remove the prior one
            distincts.get(hash).invalidate();
            distincts.remove(hash);
            redundant=false;
          }
          else redundant=true;
        }

      if (redundant==null && seqs.row(H.getProteinKey())!=null &&
          Strs.hasSubStr(H.getSequence(), seqs.row(H.getProteinKey()).keySet())) redundant=true;

      // deposite the new and distinct hit
      if (redundant!=null && redundant) E.getValue().invalidate();
      else
      {
        Tools.putKeysVal(H, distincts, hashes.get(0), hashes.get(1), mzHash);
        seqs.put(H.getProteinKey(), E.getValue().getSequence(), H);
      }
    }

    return hits;
  }

  public Ms2Hits_ consolidateByFragments(double max_miss_pct)
  {
    Table<Integer, Ms2Hit_, Set<Double>> hit_frag = HashBasedTable.create();
    if (Tools.isSet(getExactMatches())) consolidateByFragments(getExactMatches(), max_miss_pct, hit_frag);
    if (Tools.isSet(getOpenMatches( ))) consolidateByFragments(getOpenMatches( ), max_miss_pct, hit_frag);

    return purge();
  }
  public Multimap<Double, Ms2Hit_> consolidateByFragments(Multimap<Double, Ms2Hit_> hits, double max_miss_pct, Table<Integer, Ms2Hit_, Set<Double>> n_hit_frag)
  {
    if (!Tools.isSet(hits)) return hits;

    Set<Double> frags = new HashSet<>();
    for (Double scr : hits.keySet())
      for (Ms2Hit_ H : hits.get(scr))
      {
        frags.clear();
        for (AnnotatedPeak pk : H.getY().getTrack()) frags.add(pk.getIntensity());
        for (AnnotatedPeak pk : H.getB().getTrack()) frags.add(pk.getIntensity());

        // check the matching fragments to see if they're already used elsewhere
        if (!consolidateByFragment(frags, max_miss_pct, n_hit_frag))
        {
          // despoite the frags for later use
          n_hit_frag.put(frags.size(), H, new HashSet<>(frags));
        }
        else
          H.invalidate();
      }

    return hits;
  }
  private boolean consolidateByFragment(Set<Double> frags, double max_miss_pct, Table<Integer, Ms2Hit_, Set<Double>> n_hit_frag)
  {
    if (Tools.isSet(frags) && Tools.isSet(n_hit_frag))
      for (int i : n_hit_frag.rowKeySet())
      {
        int max_miss = (int )Math.round(i*max_miss_pct*0.01d);
        if (i-frags.size()>=max_miss)
          for (Ms2Hit_ E : n_hit_frag.row(i).keySet())
          {
            Sets.SetView<Double> intersect = Sets.intersection(frags, n_hit_frag.get(i, E));
            if (frags.size()==intersect.size() && (n_hit_frag.get(i, E).size()-intersect.size())<=max_miss)
            {
              // this is a duplicate
              return true;
            }
          }
      }

    return false;
  }
  public Ms2Hits_ trimCandidates(int tops)
  {
    if (!Tools.isSet(mCandidates) || mCandidates.size()<=tops) return this;

    Iterator<Map.Entry<Double, Ms2Hit_>> itr = mCandidates.entries().iterator();
    while (itr.hasNext())
    {
      if (itr.next().getValue().getRank()>tops) itr.remove();
    }

    return this;
  }
  // MUST call the 'survey' beforehand to setup the final score
  public Ms2Hits_ qvalues()
  {
    // combine the open and exact matches
    mCandidates = TreeMultimap.create(Ordering.natural().reverse(), Ordering.natural());
    if (Tools.isSet(getExactMatches()))
      for (Ms2Hit_ H : getExactMatches().values()) mCandidates.put(H.getScore(), H);

    if (Tools.isSet(getOpenMatches()))
      for (Ms2Hit_ H : getOpenMatches().values( )) mCandidates.put(H.getScore(), H);

    // count the total number of decoys
    double decoys=0; int rank=0;
    for (Double score : getCandidates().keySet())
    {
      rank++;
      for (Ms2Hit_ H : getCandidates().get(score))
        if (H.setRank(rank).isDecoy()) decoys++;
    }
    // setup the qvalue calculation
    QVals        qvals = new QVals();
    Histogram      all = new Histogram("all");
    SimpleRegression R = new SimpleRegression(true);

    int pct10 = (int )Math.min(15, Math.max(3, 0.1d * decoys));
//    System.out.println("Score\tSurvival");
    for (Ms2Hit_ F : getCandidates().values())
    {
      qvals.put(F.getScore(), F.isDecoy());
      if (F.isDecoy())
      {
        all.add(F.getScore());
        if (R.getN()<pct10)
        {
          R.addData(F.getScore(), Math.log10((double) all.getData().size() / decoys));
        }
//        System.out.println(F.getScore()+"\t"+(double )all.getData().size()/decoys);
      }
    }
//    System.out.println("Score\tSurvival");
//    Collections.sort(all.getData(), Collections.reverseOrder());
//    for (int i=0; i<all.getData().size(); i++)
//      System.out.println(all.getData().get(i)+"\t"+((double )i/(double )all.getData().size()));

    all.generate(all.getData().size()>100?25:15).assessTruncated(1);
    if (getScoreModel(Ms2Hit_.SCR_GAP)!=null && all.getCenter()!=null) getScoreModel(Ms2Hit_.SCR_GAP).setCenter(all.getCenter());

    // calculate the q-value
    qvals.model(); setScore(Ms2Hit_.SCR_THRESHOLD, qvals.getThreshold());

    // get the expect-value per ProteinProspector
    //http://prospector.ucsf.edu/prospector/html/misc/publications/2006_ASMS_1.pdf
    setScore(Ms2Hit_.SCR_EVAL_K, R.getSlope());
    setScore(Ms2Hit_.SCR_EVAL_B, R.getIntercept());

    for (Ms2Hit_ F : getCandidates().values())
    {
      F.setScore(Ms2Hit_.SCR_EVAL, -10*(R.getSlope()*F.getScore()+R.getIntercept())); // in dB format
      F.setScore(Ms2Hit_.SCR_SNR,  Math.pow(10d, (Math.log10(1d/decoys)-R.getSlope()*F.getScore()-R.getIntercept())));
    }

    return this;
  }

  public int size()
  {
    int counts = (Tools.isSet(mExactMatches)?mExactMatches.size():0) + (Tools.isSet(mOpenMatches)?mOpenMatches.size():0);
    return counts>0?counts:(getCandidates()!=null?getCandidates().size():0);
  }

  private Ms2Hits_ initFpmEntries()
  {
    mDecoyY = new Histogram("decoy y"); mDecoyB = new Histogram("decoy_b");
    // skip the duplicated scores
    for (Double score : getCtSegments().rowKeySet())
    {
      double goods = Stats.greaterEqualThan(getCtSegments().row(score).keySet(), 0d),
              bads = getCtSegments().row(score).keySet().size()-goods, base = Math.min(goods + bads, 4);
      // reduce the counts to the base
      for (int i=0; i<Math.round(bads*(base/(goods+bads))); i++) mDecoyY.add(score);
    }
    for (Double score : getNtSegments().rowKeySet())
    {
      double goods = Stats.greaterEqualThan(getNtSegments().row(score).keySet(), 0d),
              bads = getNtSegments().row(score).keySet().size()-goods, base = Math.min(goods + bads, 4);
      // reduce the counts to the base
      for (int i=0; i<Math.round(bads*(base/(goods+bads))); i++) mDecoyB.add(score);
    }

    mDecoyB.generate(25); mDecoyY.generate(25);

    if (Tools.isSet(mDecoyY.getHistogram()))
    {
      mDecoyY.assessTruncated(4);
      setScore(Ms2Hit_.SCR_DECOY_Y, Collections.max(mDecoyY.getData()));
      setScore(Ms2Hit_.SCR_DECOY_Y0, mDecoyY.getCenter());
      setScore(Ms2Hit_.SCR_DECOY_Y1, mDecoyY.getSigma());
    }
    if (Tools.isSet(mDecoyB.getHistogram()))
    {
      mDecoyB.assessTruncated(4);
      setScore(Ms2Hit_.SCR_DECOY_B, Collections.max(mDecoyB.getData()));
      setScore(Ms2Hit_.SCR_DECOY_B0,                mDecoyB.getCenter());
      setScore(Ms2Hit_.SCR_DECOY_B1,                mDecoyB.getSigma());
    }

    return this;
  }

  public Ms2Hits_ purge()
  {
    if (Tools.isSet(getExactMatches()))
    {
      Iterator<Map.Entry<Double, Ms2Hit_>> itr = getExactMatches().entries().iterator();
      while (itr.hasNext())
        if (!Strs.isSet(itr.next().getValue().getSequence())) itr.remove();
    }
    if (Tools.isSet(getOpenMatches()))
    {
      Iterator<Map.Entry<Double, Ms2Hit_>> itr = getOpenMatches().entries().iterator();
      while (itr.hasNext())
      {
        Ms2Hit_ H = itr.next().getValue();
        if (!Strs.isSet(H.getSequence())) itr.remove();
        // do not allow any open match without some matches in both ends
        else if (H.getY()==null || !Tools.isSet(H.getY().getTrack()) ||
                 H.getB()==null || !Tools.isSet(H.getB().getTrack())) itr.remove();
      }
    }
    return this;
  }
  // setup the score models in preparation for the final composite scores
  public Ms2Hits_ survey()
  {
    mScoreModels = new HashMap<>();
    mScoreModels.put(Ms2Hit_.SCR_GAP,  new ScoreModel(Ms2Hit_.SCR_GAP));
//    mScoreModels.put(Ms2Hit_.SCR_KAI,  new ScoreModel(Ms2Hit_.SCR_KAI));
//    mScoreModels.put(Ms2Hit_.SCR_MATCH,new ScoreModel(Ms2Hit_.SCR_MATCH));

    // gather the decoy data
    if (Tools.isSet(getExactMatches()))
      for (Ms2Hit_ hit : getExactMatches().values())
      {
        if (hit.isDecoy())
          for (ScoreModel m : mScoreModels.values()) m.addExactDecoy(hit.getScore(Ms2Hit_.SCR_GAP));
      }
    if (Tools.isSet(getOpenMatches()))
      for (Ms2Hit_ hit : getOpenMatches().values())
      {
        if (hit.isDecoy())
          for (ScoreModel m : mScoreModels.values()) m.addOpenDecoy(hit.getScore(Ms2Hit_.SCR_GAP));
      }

    // set the initial counts
    // setup the offsets by the peptide len
    for (ScoreModel m : mScoreModels.values())
    {
      m.setOffset(ScoreModel.eType.exact, Tools.isSet(getExactMatches()) ? getExactMatches().size() : 0);
      m.setOffset(ScoreModel.eType.open, Tools.isSet(getOpenMatches()) ? getOpenMatches().size() : 0);
      m.model(ScoreModel.eType.exact, ScoreModel.eType.open);
    }

    // create the final scores
    return adjustGapScores();
  }
  public Ms2Hits_ adjustGapScores()
  {
    if (getScoreModel(Ms2Hit_.SCR_GAP)!=null)
    {
      // no change if there wasn't enough decoys around
      double offset = getScoreModel(Ms2Hit_.SCR_GAP).getOffset(ScoreModel.eType.open)!=null?
                      getScoreModel(Ms2Hit_.SCR_GAP).getOffset(ScoreModel.eType.open):0d;
      // gather the decoy data
      if (Tools.isSet(getExactMatches()))
        for (Ms2Hit_ hit : getExactMatches().values())
          hit.setScore(Ms2Hit_.SCR_FINAL, hit.getScore(Ms2Hit_.SCR_GAP));

      if (Tools.isSet(getOpenMatches()))
        for (Ms2Hit_ hit : getOpenMatches().values())
          hit.setScore(Ms2Hit_.SCR_FINAL, hit.getScore(Ms2Hit_.SCR_GAP)-offset);
    }

    return this;
  }
  public Ms2Hits_ consolidateNegativeResidue(float[] AAs)
  {
    consolidateNegativeResidue(getExactMatches(), AAs);
    consolidateNegativeResidue(getOpenMatches(), AAs);

    purge();

    return this;
  }
  public Multimap<Double, Ms2Hit_> consolidateNegativeResidue(Multimap<Double, Ms2Hit_> hits, float[] AAs)
  {
    if (Tools.isSet(hits))
      for (Ms2Hit_ H : hits.values())
        if (Tools.isSet(H.getMods()))
          for (Integer loc : H.getMods().keySet())
            if (AAs[H.getSequence().charAt(loc-H.getLeft())]+H.getMods().get(loc)<-0.01)
            {
              H.invalidate(); break;
            }

    return hits;
  }
  // remove the duplicates
  public Ms2Hits_ consolidate(float[] AAs, OffsetPpmTolerance tol, Range<Integer> isoErr, float deci, TreeMultimap<Float, String> blocks)
  {
    // remove the invalid hits
    purge();

    HashMap<Integer, Ms2Hit_>           distincts = new HashMap<>();
    MultiTreeTable<Integer, String, Ms2Hit_> seqs = MultiTreeTable.create();

    // remove the redundant matches
    consolidate(getExactMatches(), distincts, seqs, AAs, tol, isoErr, deci, blocks);
    consolidate(getOpenMatches(),  distincts, seqs, AAs, tol, isoErr, deci, blocks);

    purge();

    return this;
  }

  public double getCenterY() { return (hasBasis(Ms2Hit_.SCR_DECOY_Y0, Ms2Hit_.SCR_DECOY_Y1)?getBasis(Ms2Hit_.SCR_DECOY_Y0):0d); }
  public double getSigmaY( ) { return (hasBasis(Ms2Hit_.SCR_DECOY_Y0, Ms2Hit_.SCR_DECOY_Y1)?getBasis(Ms2Hit_.SCR_DECOY_Y1):1d); }
  public double getCenterB() { return (hasBasis(Ms2Hit_.SCR_DECOY_B0, Ms2Hit_.SCR_DECOY_B1)?getBasis(Ms2Hit_.SCR_DECOY_B0):0d); }
  public double getSigmaB( ) { return (hasBasis(Ms2Hit_.SCR_DECOY_B0, Ms2Hit_.SCR_DECOY_B1)?getBasis(Ms2Hit_.SCR_DECOY_B1):1d); }

  public Ms2Hits_ dispose_intermediates()
  {
//    mSpectrum=null;

    // 'ladders' build from the y and b ions, respectively
    Tools.dispose(mCtSegments);
    Tools.dispose(mNtSegments);
    Tools.dispose(mExactMatches);
    Tools.dispose(mOpenMatches);
    Tools.dispose(mDecoyY, mDecoyB);
//    Tools.dispose(mCandidates);

    return this;
  }
  @Override
  public void dispose()
  {
    Tools.dispose(mPeakStats);

    mSpectrum=null;
    // 'ladders' build from the y and b ions, respectively
    Tools.dispose(mCtSegments);
    Tools.dispose(mNtSegments);
    Tools.dispose(mExactMatches);
    Tools.dispose(mOpenMatches);
    Tools.dispose(mCandidates);
  }
}
