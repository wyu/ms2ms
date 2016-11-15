package org.ms2ms.data.ms;

import com.google.common.collect.Range;
import com.google.common.collect.TreeMultimap;
import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.ms2ms.Disposable;
import org.ms2ms.algo.Isotopes;
import org.ms2ms.algo.MsStats;
import org.ms2ms.algo.Peptides;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.util.*;

/**
 * Created by yuw on 8/7/16.
 */
public class Ms2Hit implements Comparable<Ms2Hit>, Disposable
{
  public static final String SCR_KAI    = "KaiScore";
  public static final String SCR_GAP    = "GapScore";
  public static final String SCR_FINAL  = "AdjGapScore";
  public static final String SCR_MATCH  = "MatchProb";
  public static final String SCR_GLOBAL = "GlobalMatchProb";
  public static final String SCR_COMP   = "CompositeScore";
  //public static final String SCR_YB_Z   = "ntZ+ctZ";
  public static final String SCR_Y_Z    = "ctZ";
  public static final String SCR_B_Z    = "ntZ";
  public static final String SCR_THRESHOLD  = "ThresholdByQVal";
  public static final String SCR_DELTA_PKS = "incremental gains in the matched peaks";
  public static final String SCR_DELTA_MATCH = "incremental gains in the global match prob";

  public static final String SCR_DELTA  = "DeltaScore";
  public static final String SCR_OFFSET = "ScoreOffset";
  public static final String SCR_FACTOR = "ScoreFactor";
  public static final String SCR_EVAL   = "Eval";
  public static final String SCR_EVAL_B = "Eval-Intercept";
  public static final String SCR_EVAL_K = "Eval-Slope";

  public static final String SCR_SNR   = "S/N";

  public static final String SCR_DECOY_Y = "best decoy-y";
  public static final String SCR_DECOY_B = "best decoy-b";
  public static final String SCR_DECOY_Y0 = "mean of decoy-y";
  public static final String SCR_DECOY_B0 = "mean of decoy-b";
  public static final String SCR_DECOY_Y1 = "stdev of decoy-y";
  public static final String SCR_DECOY_B1 = "stdev of decoy-b";

  public static final String TAG_EXACT  = "E";
  public static final String TAG_OPEN   = "O";

  private FpmEntry mY, mB;
  private Integer mProteinKey;
  private int mCharge, mLeft, mRight, // 0-based index of the first and last residue
              mRank, mIsotopeError=0/*, mPrecursorCharge=0*/;
  private Peak mCalc=null;
  private double mDeltaM;
  private String mSequence, mPrev, mNext, mTag;
  private TreeMap<Integer, Double> mMods;
  private Map<String, Double> mScores = new HashMap<>();

  public static class ScoreDesendComparator implements Comparator<Ms2Hit>
  {
    String score;
    public ScoreDesendComparator(String s) { score=s; }
    public int compare(Ms2Hit o1, Ms2Hit o2)
    {
      return o1 != null && o2 != null ? Double.compare(o2.getScore(score), o1.getScore(score)) : 0;
    }
  }

  public Ms2Hit()
  {
    super();
    mY = new FpmEntry(); mB = new FpmEntry();
  }
  public Ms2Hit(String sequence)
  {
    super(); mSequence=sequence;
  }
  public Ms2Hit(int protein, FpmEntry y, FpmEntry b, int left, int right, int z)
  {
    super();
    mProteinKey=protein; mY=y; mB=b; setLocation(left, right); mCharge=z;
  }

  public boolean  isDecoy()       { return mProteinKey<0; }
  public Integer  getProteinKey() { return mProteinKey; }
  public int      getLeft()       { return mLeft; }
  public int      getRight()      { return mRight; }
  public int      getRank()       { return mRank; }
  public int      getIsotopeError() { return mIsotopeError; }
  public int      getCharge()     { return mCharge; }
  public Double   getScore(String s) { return mScores.get(s); }
  public double   getDelta()      { return mDeltaM; }
  public double   getCalcMH()     { return mCalc!=null?mCalc.getMz():0; }
//  public Double   getScoreOffset(){ return mScores!=null&&mScores.get(SCR_OFFSET)!=null?mScores.get(SCR_OFFSET):0d; }
  public Double   getScore()      { return mScores.get(SCR_FINAL); } // the official score!!
//  public Double   getKaiScore()   { return mScores.get(SCR_KAI); }
//  public Double   getFactor()     { return mScores.containsKey(SCR_FACTOR)?mScores.get(SCR_FACTOR):1d; }
//  public Double   getDeltaScore() { return mScores.get(SCR_DELTA); }
//  public Double   getMatchProb()  { return mScores.get(SCR_MATCH); }
  public FpmEntry getY()          { return mY; }
  public FpmEntry getB()          { return mB; }
//  public Double   getGapScore()   { return ((mY!=null?mY.getGapScore():0)+(mB!=null?mB.getGapScore():0))/getFactor()+getScoreOffset(); }
  public Double   getGapScore()   { return  (mY!=null?mY.getGapScore():0)+(mB!=null?mB.getGapScore():0); }
  public String   getPeptide()    { return ((mPrev!=null?mPrev.charAt(mPrev.length()-1):"-")+"."+mSequence+"."+(mNext!=null?mNext.charAt(0):"-")); }
  public String   getSequence()   { return mSequence; }
  public String   getPrev()       { return mPrev; }
  public String   getNext()       { return mNext; }
  public String   getTag()        { return mTag; }

  public double   getModMass()    { return Tools.isSet(mMods)?Stats.sum(mMods.values()):0d; }
//  public double   getScore(Map<String, Double> basis, double w)
//  {
//    // a composite score of several components
//    double gap = (getGapScore( )-basis.get(Ms2Hits.CTR_GAP  ))/(basis.get(Ms2Hits.SIG_GAP  )!=null?basis.get(Ms2Hits.SIG_GAP  ):1d),
//         match = (getMatchProb()-basis.get(Ms2Hits.CTR_MATCH))/(basis.get(Ms2Hits.SIG_MATCH)!=null?basis.get(Ms2Hits.SIG_MATCH):1d);
//
//    return gap+w*match;
//  }
//  public double score(Map<String, ScoreModel> models, ScoreModel.eType type)
//  {
//    // a composite score of several components
//    double scr=0d, ws=0d;
//    scr+=models.get(SCR_GAP  ).score(getScore(SCR_GAP  ),type); ws+=models.get(SCR_GAP  ).getWeight();
//    scr+=models.get(SCR_MATCH).score(getScore(SCR_MATCH),type); ws+=models.get(SCR_MATCH).getWeight();
//    scr+=models.get(SCR_KAI  ).score(getScore(SCR_KAI  ),type); ws+=models.get(SCR_KAI  ).getWeight();
//
//    // normalize by the peptide length, scaled to 100d
//    setScore(SCR_COMP, 10d*scr/(ws));
//    return getComposite();
//  }
//  @Deprecated
//  public double updateScore(Map<String, Double> basis, double w_match, double w_kai)
//  {
//    // a composite score of several components
//    double scr=0d, ws=0d;
//    if (Tools.hasKeys(basis, Ms2Hits.CTR_GAP,   Ms2Hits.SIG_GAP  )) { scr +=        (getGapScore( )-basis.get(Ms2Hits.CTR_GAP  ))/basis.get(Ms2Hits.SIG_GAP);   ws+=1d; };
//    if (Tools.hasKeys(basis, Ms2Hits.CTR_MATCH, Ms2Hits.SIG_MATCH)) { scr +=w_match*(getMatchProb()-basis.get(Ms2Hits.CTR_MATCH))/basis.get(Ms2Hits.SIG_MATCH); ws+=w_match; };
//    if (Tools.hasKeys(basis, Ms2Hits.CTR_KAI,   Ms2Hits.SIG_KAI  )) { scr +=w_kai  *(getKaiScore()-basis.get(Ms2Hits.CTR_KAI   ))/basis.get(Ms2Hits.SIG_KAI);   ws+=w_kai; };
//
//    setScore(SCR_COMP, scr/ws);
//    return getComposite();
//  }
//  public TreeMap<Integer, Double> getMods() { return mMods; }

  public Map<Integer, Double> getMods() { return mMods; }

  public Map<Integer, Double> getMod0()
  {
    SortedMap<Integer, Double> out = new TreeMap<>();
    if (Tools.isSet(mMods) && Collections.max(mMods.keySet())>=getLeft())
    {
      for (Map.Entry<Integer, Double> E : mMods.entrySet())
        if (E.getKey()>=getLeft() && E.getKey()<=getRight())
          out.put(E.getKey()-getLeft(), E.getValue());
    }

    return out;
  }
  public void invalidate()
  {
    mSequence=null; mY=null; mB=null;
  }
//  public Ms2Hit   setPrecursorCharge(int s)
//  {
//    mCharge=s
//    if (mCalc==null) mCalc = new Peak();
//    mCalc.setMzAndCharge(mCalc.getMz(), s); return this;
//  }
//  public Ms2Hit   setLeft(          int s) { mLeft =s; return this; }
//  public Ms2Hit   setRight(         int s) { mRight=s; return this; }
  public Ms2Hit   setLocation(int left,int right)
  {
//    if (left>=0 && right>=0 && left>right)
//      System.out.print("");

    mRight=right; mLeft=left; return this;
  }
  public Ms2Hit   setCharge(        int s) { mCharge=s; return this; }
  public Ms2Hit   setRank(          int s) { mRank=s; return this; }
  //public Ms2Hit   setPeptide(String s) { mPeptide=s; return this; }
  public Ms2Hit   setScore(String k, Double s) { mScores.put(k,s); return this; }
//  public Ms2Hit   setGapScore()            { mScores.put(SCR_GAP, (mY!=null?mY.getGapScore():0)+(mB!=null?mB.getGapScore():0)); return this; }
  public Ms2Hit   setGapScore()
  {
    mScores.put(SCR_GAP, getGapScore()); return this;
  }
//  public Ms2Hit   setEVal(       double s) { mScores.put(SCR_EVAL,   s); return this; }
  public Ms2Hit   setKaiScore(   double s) { mScores.put(SCR_KAI,-1d*s); return this; }
//  public Ms2Hit   setDeltaScore( double s) { mScores.put(SCR_DELTA,  s); return this; }
//  public Ms2Hit   setScoreOffset(double s) { mScores.put(SCR_OFFSET, s); return this; }
//  public Ms2Hit   setFactor(     double s) { mScores.put(SCR_FACTOR, s); return this; }
//  public Ms2Hit   setMatchProb(  double s) { mScores.put(SCR_MATCH,  s); return this; }
//  public Ms2Hit   setTag(        String s) { mTag=s; return this; }
  public Ms2Hit   extTag(        String s) { mTag=Strs.extend(mTag,s,"."); return this; }

  public Ms2Hit   setDelta(      double s) { mDeltaM=s; return this; }
  public Ms2Hit   setMH(double calc, double delta)
  {
    if (mCalc==null) mCalc = new Peak();

    mCalc.setMzAndCharge(calc, mCalc.getCharge()); mDeltaM=delta;
    return this;
  }
  public Ms2Hit   increMH(double m)
  {
    if (mCalc==null) mCalc = new Peak();
    mCalc.setMzAndCharge(mCalc.getMz()+m,mCalc.getCharge()); mDeltaM-=m; return this;
  }
  public Ms2Hit setNext(String s) { mNext=s; return this; }
  public Ms2Hit setPrev(String s) { mPrev=s; return this; }

  public  Ms2Hit setSequence(String s) { mSequence=s; return this; }
  public  Ms2Hit setPeptide(String sequence) { return sequence!=null?setPeptide(sequence.toCharArray(), 0, sequence.length()-1):this; }
  public  Ms2Hit setPeptide(char[] sequence) { return setPeptide(sequence, getLeft(), getRight()); }
  public Ms2Hit setPeptide(char[] sequence, int left, int right)
  {
    if (sequence!=null && left>=0 && right<sequence.length && right>left)
    {
      setLocation(left, right);
      mPrev    =(getLeft()>0?Strs.toString(sequence, Math.max(0, getLeft() - 5), getLeft()):"-");
      mSequence=Strs.toString(sequence, getLeft(), getRight()+1); // 0-based index
      mNext    =(getRight()<sequence.length-1?Strs.toString(sequence, getRight()+1, Math.min(sequence.length, getRight()+6)):"-");
    }
    return this;
  }
  public Ms2Hit clearMods() { Tools.dispose(mMods); return this; }
  public Ms2Hit purgeMods()
  {
    // remove the mod if it's too small
    if (Tools.isSet(mMods))
    {
      Iterator<Map.Entry<Integer, Double>> itr = mMods.entrySet().iterator();
      while (itr.hasNext())
      {
        Double m = itr.next().getValue();
        if (Math.abs(m)<0.1) { setDelta(getDelta()+m); itr.remove(); }
      }
      if (mMods.size()==0) mMods=null;
    }
    return this;
  }
  public boolean hasMod(int loc) { return mMods!=null&&mMods.containsKey(loc); }
  public Ms2Hit setMod(int loc, double mod)
  {
    if (mMods==null) mMods = new TreeMap<>(); else mMods.clear();
    mMods.put(loc, mod); return this;
  }
  public Ms2Hit addMod(int loc, double mod)
  {
    if (mMods==null) mMods = new TreeMap<>();

    if (mMods.containsKey(loc) && !mMods.get(loc).equals(mod))
    {
//      System.out.println("        Conflicting mod from N/C terminus!");
    }
    // enforcing the terminal positions, Nov 4, 2016
    else if (loc>=getLeft() && loc<=getRight()) mMods.put(loc, mod);

    return this;
  }
  public Ms2Hit increMod(int loc, double mod)
  {
    if (mMods==null) mMods = new TreeMap<>();

    // increment the mod at the same location
    if (mMods.containsKey(loc))     mod += mMods.get(loc);
    if (loc>=getLeft() && loc<=getRight()) mMods.put(loc, mod);

    return this;
  }
  public Ms2Hit   setModLocation(List<ModLocation> N, List<ModLocation> C)
  {
    if      ( Tools.isSet(N) && !Tools.isSet(C)) addModLocation(N);
    else if (!Tools.isSet(N) &&  Tools.isSet(C)) addModLocation(C);
    else if ( Tools.isSet(N) &&  Tools.isSet(C))
    {
      // better to be consistent with each other
      addModLocation(N).addModLocation(C);
    }

    return this;
  }
  private Ms2Hit addModLocation(Collection<ModLocation> mods)
  {
    if (Tools.isSet(mods))
    {
      if (mMods==null) mMods = new TreeMap<>();
      for (ModLocation mod : mods)
      {
        addMod(mod.locations, mod.mods);
//        if (mMods.containsKey(mod.locations) && !mMods.get(mod.locations).equals(mod.mods))
//          throw new RuntimeException("Conflicting mod from N/C terminus!");
//        else mMods.put(mod.locations, mod.mods);
      }
    }

    return this;
  }
  public Ms2Hit setFpmEntries(FpmEntry y, FpmEntry b) { mY=null; mY=y; mB=null; mB=b; return this; }

//  public Ms2Hit setProb()
//  {
//    mProb=(getY()!=null?getY().getProb():0d)+(getB()!=null?getB().getProb():0d);
//    return this;
//  }
//  public Ms2Hit calcScore(double precision_ppm)
//  {
//    NormalDistribution norm = new NormalDistributionImpl(0, precision_ppm);
//    // TODO to be completed
//
//    return this;
//  }
  public Ms2Hit shallow_copy()
  {
    Ms2Hit clone = new Ms2Hit(getProteinKey(), mY.shallow_copy(), mB.shallow_copy(), mLeft, mRight, mCharge);
    clone.setMH(getCalcMH(), mDeltaM).setPeptide(mSequence);

    return clone;
  }
  @Override
  public String toString()
  {
    String mods="";
    if (Tools.isSet(mMods))
      for (Map.Entry<Integer, Double> E : mMods.entrySet())
        mods = Strs.extend(mods, Tools.d2s(E.getValue(), 3)+"@"+E.getKey(),";");

    if (Strs.isSet(mods)) mods = "("+mods+")";

    return getProteinKey()+":"+getLeft()+"-"+getRight()+",m/z"+ Tools.d2s(getCalcMH(), 5)+"$"+
        (mB!=null?mB.getTrack().size():"*")+"->"+(mY!=null?mY.getTrack().size():"*")+(Strs.isSet(getSequence())?("="+getPeptide()):"")+mods+"^"+
         MsStats.asDeviation(mDeltaM, getCalcMH(), 999) + "->" + Tools.d2s(getGapScore(), 2);
  }

  @Override
  public int compareTo(Ms2Hit o)
  {
    int c = Stats.compareTo(mProteinKey, o.getProteinKey());

    if (c==0) c = Stats.compareTo(mLeft, o.getLeft());
    if (c==0) c = Stats.compareTo(mRight, o.getRight());
    if (c==0)
    {
      // check the mods
      int mod0=(mMods!=null?mMods.size():0), mod1=(o.mMods!=null?o.mMods.size():0);
      c = Integer.compare(mod0, mod1);
      // don;t really care who is ahead now
      if (c==0 && mod0>0) c = Integer.compare(mMods.hashCode(), o.mMods.hashCode());
    }

    return c;
  }
// check if the MH delta is consistent with isotope error following MSGF+
//  [-ti IsotopeErrorRange] (Range of allowed isotope peak errors, Default:0,1)
//  Takes into account of the error introduced by chooosing a non-monoisotopic peak for fragmentation.
//  E.g. "-t 20ppm -ti -1,2" tests abs(exp-calc-n*1.00335Da)<20ppm for n=-1, 0, 1, 2.
  public Ms2Hit isIsotopeErr(Range<Integer> range, Tolerance tol)
  {
    if (Tools.isSet(range))
      for (int i=range.lowerEndpoint(); i<=range.upperEndpoint(); i++)
//        if (tol.withinTolerance(getDelta()+getCalcMH(), i*Isotopes.DELTA_C13+getCalcMH()))
        if (tol.withinTolerance(getCalcMH(), i*Isotopes.DELTA_C13+getCalcMH()+getDelta()))
        {
          mIsotopeError=i; mDeltaM-=(i*Isotopes.DELTA_C13);
          break;
        }

    return this;
  }
  public Integer hashcodeByYBZ()
  {
    int hash=0;
    if (getY()!=null) hash+=getY().hashcodeByTrack();
    if (getB()!=null) hash+=getB().hashcodeByTrack();

    return hash+getCharge();
  }
  // by the observed m/z only, not position
  public Integer hashcodeByYBmz()
  {
    int hash=0;
    if (getY()!=null) hash+=getY().hashcodeByTrackMz();
    if (getB()!=null) hash+=getB().hashcodeByTrackMz();

    return hash;
  }
  //  public Integer hashcodeByIntervals(float[] AAs)
//  {
//    // straight up hash, not worry about the mods, etc
//    if (Strs.isSet(getSequence()))
//    {
//      float[] intervals = new float[getSequence().length()];
//      for (int i=0; i<getSequence().length(); i++) intervals[i]=AAs[getSequence().charAt(i)];
//
//      int hash=0;
//      for (int i=0; i<intervals.length; i++) hash+=(i+1)*intervals[i]/tol;
//      return hash+getCharge();
//    }
//
//    return getCharge();
//  }
//  public List<Integer> hashcodeByIntervals(float[] AAs, OffsetPpmTolerance tolerance, Range<Integer> isotopeErr, float deci)
//  {
//    if (getCalcMH()==0) mCalc=new Peak(Peptides.calcMH(getSequence().toCharArray(), 0, getSequence().length()-1, AAs), 0d);
//
//    float tol = (float )(tolerance.getMax(getCalcMH())-tolerance.getMin(getCalcMH()));
//    if (Strs.isSet(getSequence()))
//    {
//      List<Float> intervals = new ArrayList<>();
//      for (int i=0; i<getSequence().length(); i++) intervals.add(AAs[getSequence().charAt(i)]);
//
//      // let's consider the case of false-extension
//      List<Integer> putatives = null;
//      Integer  mod = Tools.isSet(mMods)?Collections.max(mMods.keySet()):null;
//      double delta = Tools.isSet(mMods)?mMods.get(mod):0;
//      if (     delta<-56)
//      {
//        putatives = Peptides.seekRemoval(getSequence(), true,  getCalcMH(), delta, tolerance, isotopeErr, AAs);
//        if (Tools.isSet(putatives))
//        {
//          mMods.remove(mod);
//          // TODO to be tested
//          for (Integer putative : putatives)
//            intervals.remove(putative.intValue());
//        }
//      }
//      else if (delta> 56)
//      {
//        putatives = Peptides.seekRemoval(getNext(), false, getCalcMH(), delta, tolerance, isotopeErr, AAs);
//        if (Tools.isSet(putatives))
//        {
//          mMods.remove(mod);
//          for (Integer putative : putatives) intervals.add(AAs[getNext().charAt(putative)]);
//        }
//      }
//
//      // need the base one without mod
//      int hash0=0;
//      for (int i=0; i<intervals.size(); i++) hash0+=(i+1)*intervals.get(i)/deci;
//
//      List<List<Float>> Intervals = new ArrayList<>(), Additionals = new ArrayList<>();
//
//      // now the increment due to site-specific mods
//      if (Tools.isSet(getMod0()))
//      {
//        // check if the mod==AA
//        for (Integer m : getMod0().keySet())
//        {
//          Float AA = Tools.findClosest(AAs, getMod0().get(m).floatValue(), deci);
//          if (AA!=null)
//          {
//            List<Float> cloned = new ArrayList<>(intervals);
//            cloned.add(m+1, AA); Additionals.add(cloned);
//            List<Float> clone2 = new ArrayList<>(intervals);
//            clone2.add(m,   AA); Additionals.add(clone2);
//          }
//        }
//        for (Integer m : getMod0().keySet())
//        {
//          intervals.set(m, intervals.get(m) + getMod0().get(m).floatValue());
//          // replace it with the accurate mass if matches to an AA
//          Float AA = Tools.findClosest(AAs, intervals.get(m), deci);
//          if (AA != null) intervals.set(m, AA);
//        }
//      }
//      // trim the residue if the mass is zero within the tolerance
//      Set<Integer> removed = new HashSet<>();
//      for (int i=0; i<intervals.size(); i++)
//      {
//        if      (                        Math.abs(intervals.get(i))               <=tol)   removed.add(i);
////        if      (i>0 &&                  Math.abs(intervals[i-1])             <=tol)   removed.add(i-1); // in case the localization is not perfect
////        if      (i<intervals.length-1 && Math.abs(intervals[i+1])             <=tol)   removed.add(i+1); // in case the localization is not perfect
//        else if (i>0 &&                  Math.abs(intervals.get(i-1)+intervals.get(i))<=tol) { removed.add(i-1); removed.add(i); }
//        else if (i<intervals.size()-1 && Math.abs(intervals.get(i)  +intervals.get(i+1))<=tol) { removed.add(i);   removed.add(i+1); }
//      }
//
//      Intervals.add(intervals); Intervals.addAll(Additionals);
//
//      // enumerate the situation where neighboring residues can be combined to one
//      if (!Tools.isSet(removed) && Tools.isSet(mMods) && mMods.size()==1)
//      {
//        int pos = Tools.front(mMods.keySet())-getLeft();
//        if (pos>0 && pos<intervals.size())
//        {
//          // with residue prior
//          float m = intervals.get(pos)+intervals.get(pos-1);
//          Collection<Float> found = Tools.find(AAs, m - deci, m + deci);
//          if (Tools.isSet(found))
//          {
//            // remove the old residues
//            List<Float> cloned = new ArrayList<>(intervals);
//            cloned.remove(pos-1); cloned.remove(pos - 1);
//            for (Float A : found)
//            {
//              List<Float> c = new ArrayList<>(cloned);
//              c.add(pos-1, A);
//              Intervals.add(c);
//            }
//          }
//        }
//        if (pos>=0 && pos<intervals.size()-1)
//        {
//          // with residue after
//          float m = intervals.get(pos)+intervals.get(pos+1);
//          Collection<Float> found = Tools.find(AAs, m - deci, m + deci);
//          if (Tools.isSet(found))
//          {
//            // remove the old residues
//            List<Float> cloned = new ArrayList<>(intervals);
//            cloned.remove(pos); cloned.remove(pos);
//            for (Float A : found)
//            {
//              List<Float> c = new ArrayList<>(cloned);
//              c.add(pos, A);
//              Intervals.add(c);
//            }
//          }
//        }
//        // now consider multiple residue for one K
//        // TTTIGA(+9.1134)NY
//        // TTTIGK
//        if (pos>=0 && pos<intervals.size()-1)
//        {
//          float m=0; // run it up to the C-t end
//          for (int k=pos; k<intervals.size(); k++) m+=intervals.get(k);
//          if (Math.abs(m-AAs['K'])<=deci)
//          {
//            // remove the old residues
//            List<Float> cloned = new ArrayList<>(intervals);
//            for (int k=pos; k<intervals.size(); k++) cloned.remove(pos);
//            cloned.add(pos, AAs['K']);
//            // another possibility
//            Intervals.add(cloned);
//          }
//        }
//      }
//      List<Integer> hashes = new ArrayList<>();
//
//      int hash=0, k=0;
//      for (int i=0; i<intervals.size(); i++)
//      {
//        // no need to remove it, only out of hash
//        if (!Tools.isSet(removed) || !removed.contains(i)) { hash+=Math.round((k+1)*intervals.get(i)/deci); k++; }
//      }
//      hashes.add(hash+getCharge()); hashes.add(hash0+getCharge());
//
//      if (Intervals.size()>1)
//      for (int j=1; j<Intervals.size(); j++)
//      {
//        hash=k=0;
//        for (int i=0; i<Intervals.get(j).size(); i++)
//        {
//          // no need to remove it, only out of hash
//          hash+=Math.round((i+1)*Intervals.get(j).get(i)/deci);
//        }
//        hashes.add(hash+getCharge());
//      }
//
//      return hashes;
//    }
//
//    return null;
//  }
  // "blocks" is a dictionary of AA combination to be considered as building blocks
  public List<Integer> hashcodeByIntervals(float[] AAs, OffsetPpmTolerance tolerance, Range<Integer> isotopeErr, float deci, TreeMultimap<Float, String> blocks)
  {
    // everything is built on the peptide sequence
    if (!Strs.isSet(getSequence())) return null;

    // populate the calculate MH if necessary
    if (getCalcMH()==0) mCalc=new Peak(Peptides.calcMH(getSequence().toCharArray(), 0, getSequence().length()-1, AAs), 0d);

    // setup the basic mass interval with just the AA sequence, not including the mods
    List<Integer>  hashes = new ArrayList<>();
    List<Float> intervals = new ArrayList<>(); int hash0=0; // hash without mods
    for (int i=0; i<getSequence().length(); i++)
    {
      intervals.add(AAs[getSequence().charAt(i)]);
      hash0+=Math.round((i+1)*intervals.get(i)/deci);
    }

    // quit if we don't have any mod to consider
    if (!Tools.isSet(mMods))
    {
      hashes.add(hash0+getCharge()); return hashes;
    }

    // let's consider the case of false-extension
    float    tol = (float )(tolerance.getMax(getCalcMH())-tolerance.getMin(getCalcMH()));
    Integer  mod = Collections.max(mMods.keySet());
    double delta = mMods.get(mod);

    List<Integer> putatives;
    if (delta<-56 || delta>56)
    {
      putatives = Peptides.seekRemoval(delta<-56?getSequence():getNext(), delta<-56,getCalcMH(),delta,tolerance,isotopeErr,AAs);
      if (Tools.isSet(putatives))
      {
        mMods.remove(mod);
        for (Integer putative : putatives)
          if (delta<-56) intervals.remove(putative.intValue()); else intervals.add(AAs[getNext().charAt(putative)]);
      }
    }

    List<List<Float>> Additionals = new ArrayList<>();

    // now the increment due to site-specific mods to check if the mod==AA
    for (Integer m : getMod0().keySet())
    {
      Float AA = Tools.findClosest(AAs, getMod0().get(m).floatValue(), deci);
      if (AA!=null)
      {
        List<Float> cloned = new ArrayList<>(intervals);
        cloned.add(m+1, AA); Additionals.add(cloned);
        List<Float> clone2 = new ArrayList<>(intervals);
        clone2.add(m,   AA); Additionals.add(clone2);
      }
      // increment the residue mass by the mod
      intervals.set(m, intervals.get(m) + getMod0().get(m).floatValue());
      // replace it with the accurate mass if matches to an AA
      AA = Tools.findClosest(AAs, intervals.get(m), deci);
      if (AA != null) intervals.set(m, AA);
    }

    // trim the residue if the mass is zero within the tolerance
    Set<Integer> removed = new HashSet<>();
    for (int i=0; i<intervals.size(); i++)
    {
      if      (                        Math.abs(intervals.get(i))                     <=tol)   removed.add(i);
      else if (i>0 &&                  Math.abs(intervals.get(i-1)+intervals.get(i))  <=tol) { removed.add(i-1); removed.add(i); }
      else if (i<intervals.size()-1 && Math.abs(intervals.get(i)  +intervals.get(i+1))<=tol) { removed.add(i);   removed.add(i+1); }
    }

    if (!Tools.isSet(removed) && mMods.size()==1)
    {
      // enumerate the situation where neighboring residues can be combined to one
      int pos = Tools.front(mMods.keySet())-getLeft();
      if (pos>0 && pos<intervals.size())
      {
        // with residue prior
        float m = intervals.get(pos)+intervals.get(pos-1);
        Collection<Float> found = Tools.find(AAs, m - deci, m + deci);
        if (Tools.isSet(found))
        {
          // remove the old residues
          List<Float> cloned = new ArrayList<>(intervals);
          cloned.remove(pos-1); cloned.remove(pos - 1);
          for (Float A : found)
          {
            List<Float> c = new ArrayList<>(cloned);
            c.add(pos-1, A);
            Additionals.add(c);
          }
        }
      }
      if (pos>=0 && pos<intervals.size()-1)
      {
        // with residue after
        float m = intervals.get(pos)+intervals.get(pos+1);
        Collection<Float> found = Tools.find(AAs, m - deci, m + deci);
        if (Tools.isSet(found))
        {
          // remove the old residues
          List<Float> cloned = new ArrayList<>(intervals);
          cloned.remove(pos); cloned.remove(pos);
          for (Float A : found)
          {
            List<Float> c = new ArrayList<>(cloned);
            c.add(pos, A);
            Additionals.add(c);
          }
        }

        // now consider multiple residue for one K
        // TTTIGA(+9.1134)NY
        // TTTIGK
        m=0; // run it up to the C-t end
        for (int k=pos; k<intervals.size(); k++) m+=intervals.get(k);
        if (Math.abs(m-AAs['K'])<=deci)
        {
          // remove the old residues
          List<Float> cloned = new ArrayList<>(intervals);
          for (int k=pos; k<intervals.size(); k++) cloned.remove(pos);
          cloned.add(pos, AAs['K']);
          // another possibility
          Additionals.add(cloned);
        }
      }
      if (pos>=0 && pos<intervals.size())
      {
        // another scenario where the mod'd residue can be a 2-residues combination
        Map<Float, Collection<String>> slice = blocks.asMap().subMap(intervals.get(pos)-deci, intervals.get(pos)+deci);
        if (Tools.isSet(slice))
          for (Float mm : slice.keySet())
            for (String AA2 : slice.get(mm))
            {
              List<Float> cloned = new ArrayList<>(intervals);
              cloned.remove(pos);
              for (Character c : AA2.toCharArray()) cloned.add(pos, AAs[c]);
              // another possibility
              Additionals.add(cloned);
            }
      }
    }
    // TODO something more outlandish: if the first or last two residues are not defined by the b/y ions, check the other permutations


    // output the hashes from each possibilities
    int hash=0, k=0;
    for (int i=0; i<intervals.size(); i++)
    {
      // no need to remove it, only out of hash
      if (!Tools.isSet(removed) || !removed.contains(i)) { hash+=Math.round((k+1)*intervals.get(i)/deci); k++; }
    }
    hashes.add(hash+getCharge()); hashes.add(hash0+getCharge());

    // including the alternatives
    if (Tools.isSet(Additionals))
      for (int j=0; j<Additionals.size(); j++)
      {
        hash=k=0;
        for (int i=0; i<Additionals.get(j).size(); i++)
        {
          // no need to remove it, only out of hash
          hash+=Math.round((i+1)*Additionals.get(j).get(i)/deci);
        }
        hashes.add(hash+getCharge());
      }

    return hashes;
  }
  @Override
  public void dispose()
  {
    Tools.dispose(mY, mB);
    mSequence=mPrev=mNext=null;
    Tools.dispose(mMods);
  }
}
