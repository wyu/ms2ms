package org.ms2ms.data.ms;

import com.google.common.collect.Range;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.Disposable;
import org.ms2ms.algo.MsStats;
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
  public static final String SCR_DELTA  = "DeltaScore";
  public static final String SCR_OFFSET = "ScoreOffset";
  public static final String SCR_MATCH  = "MatchProb";
  public static final String SCR_EVAL   = "Eval";
  public static final String SCR_COMP   = "CompositeScore";

  private FpmEntry mY, mB;
  private Long mProteinKey;
  private int mLeft, mRight, // 0-based index of the first and last residue
              mRank, mIsotopeError=0, mPrecursorCharge=0;
  private double mCalcMH, mDeltaM;
  //private double mKaiScore, mDeltaScore, mEval, mScoreOffset=0, mMatchProb, mGapQval, mMatchQval;
  private String mSequence, mPrev, mNext;
  private TreeMap<Integer, Double> mMods;
  private Map<String, Double> mScores = new HashMap<>();

  public Ms2Hit()
  {
    super();
    mY = new FpmEntry(); mB = new FpmEntry();
  }
  public Ms2Hit(Long protein, FpmEntry y, FpmEntry b, int left, int right)
  {
    super();
    mProteinKey=protein; mY=y; mB=b; setLocation(left, right);
  }

  public boolean  isDecoy()       { return mProteinKey!=null && mProteinKey<0; }
  public Long     getProteinKey() { return mProteinKey; }
  public int      getLeft()       { return mLeft; }
  public int      getRight()      { return mRight; }
  public int      getRank()       { return mRank; }
  public double   getScore(String s) { return mScores.get(s); }
  public double   getDelta()      { return mDeltaM; }
  public double   getCalcMH()     { return mCalcMH; }
  public Double   getEVal()       { return mScores.get(SCR_EVAL); }
  public Double   getScoreOffset(){ return mScores!=null&&mScores.get(SCR_OFFSET)!=null?mScores.get(SCR_OFFSET):0d; }
//  public double   getGapQval()    { return mGapQval; }
//  public double   getMatchQval()  { return mMatchQval; }
  public Double   getKaiScore()   { return mScores.get(SCR_KAI); }
  public Double   getDeltaScore() { return mScores.get(SCR_DELTA); }
  public Double   getMatchProb()  { return mScores.get(SCR_MATCH); }
  public FpmEntry getY()          { return mY; }
  public FpmEntry getB()          { return mB; }
//  public int      getMotifs()     { return(mY!=null?mY.getMotifs(  ):0)+(mB!=null?mB.getMotifs(  ):0); }
  public Double   getGapScore()   { return(mY!=null?mY.getGapScore():0)+(mB!=null?mB.getGapScore():0)+getScoreOffset(); }
//  public Double   getGapScore()   { return mScores.get(SCR_GAP)+getScoreOffset(); }
  public String   getPeptide()    { return (mPrev+"."+mSequence+"."+mNext); }
  public String   getSequence()   { return mSequence; }
  public double   getModMass()    { return mMods!=null? Stats.sum(mMods.values()):0d; }
  public double   getScore(Map<String, Double> basis)
  {
    // a composite score of several components
    double gap = (getGapScore( )-basis.get(Ms2Hits.CTR_GAP  ))/basis.get(Ms2Hits.SIG_GAP),
         match = (getMatchProb()-basis.get(Ms2Hits.CTR_MATCH))/basis.get(Ms2Hits.SIG_MATCH);

    return gap+0.25*match;
  }
//  public TreeMap<Integer, Double> getMods() { return mMods; }

  public Map<Integer, Double> getMod0()
  {
    if (mMods==null) return null;

    SortedMap<Integer, Double> out = new TreeMap<>();
    for (Map.Entry<Integer, Double> E : mMods.entrySet())
      out.put(E.getKey()-getLeft(), E.getValue());

    return out;
  }

  public Ms2Hit   setPrecursorCharge(int s) { mPrecursorCharge=s; return this; }
//  public Ms2Hit   setLeft(          int s) { mLeft =s; return this; }
//  public Ms2Hit   setRight(         int s) { mRight=s; return this; }
  public Ms2Hit   setLocation(int left,int right)
  {
//    if (left>=0 && right>=0 && left>right)
//      System.out.print("");

    mRight=right; mLeft=left; return this;
  }
  public Ms2Hit   setRank(          int s) { mRank=s; return this; }
  //public Ms2Hit   setPeptide(String s) { mPeptide=s; return this; }
  public Ms2Hit   setScore(String k, Double s) { mScores.put(k,s); return this; }
//  public Ms2Hit   setGapScore()            { mScores.put(SCR_GAP, (mY!=null?mY.getGapScore():0)+(mB!=null?mB.getGapScore():0)); return this; }
  public Ms2Hit   setGapScore()            { mScores.put(SCR_GAP, getGapScore()); return this; }
  public Ms2Hit   setEVal(       double s) { mScores.put(SCR_EVAL,   s); return this; }
  public Ms2Hit   setKaiScore(   double s) { mScores.put(SCR_KAI,    s); return this; }
  public Ms2Hit   setDeltaScore( double s) { mScores.put(SCR_DELTA,  s); return this; }
  public Ms2Hit   setScoreOffset(double s) { mScores.put(SCR_OFFSET, s); return this; }
  public Ms2Hit   setMatchProb(  double s) { mScores.put(SCR_MATCH,  s); return this; }
  public Ms2Hit   setMH(double calc, double delta)
  {
    mCalcMH=calc; mDeltaM=delta;
    return this;
  }
  public  Ms2Hit setPeptide(String sequence) { return sequence!=null?setPeptide(sequence.toCharArray(), 0, sequence.length()-1):this; }
  public  Ms2Hit setPeptide(char[] sequence) { return setPeptide(sequence, getLeft(), getRight()); }
  private Ms2Hit setPeptide(char[] sequence, int left, int right)
  {
    setLocation(left, right);
    mPrev    =(getLeft()>0?sequence[getLeft()-1]:'-')+"";
    mSequence=Strs.toString(sequence, getLeft(), getRight()+1); // 0-based index
    mNext    =(getRight()<sequence.length-1?sequence[getRight()+1]:'-')+"";

//    if (Tools.equals("AIAVRS", mSequence))
//      System.out.print("");
    return this;
  }
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
      System.out.println("        Conflicting mod from N/C terminus!");
    }
    else mMods.put(loc, mod);

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
  public Ms2Hit calcScore(double precision_ppm)
  {
    NormalDistribution norm = new NormalDistributionImpl(0, precision_ppm);
    // TODO to be completed

    return this;
  }
  public Ms2Hit shallow_copy()
  {
    Ms2Hit clone = new Ms2Hit(getProteinKey(), mY.shallow_copy(), mB.shallow_copy(), mLeft, mRight);
    clone.setMH(mCalcMH, mDeltaM).setPeptide(mSequence);

    return clone;
  }
  @Override
  public String toString()
  {
    return getProteinKey()+":"+getLeft()+"-"+getRight()+",m/z"+ Tools.d2s(mCalcMH, 5)+"$"+
        (mB!=null?mB.getTrack().size():"*")+"->"+(mY!=null?mY.getTrack().size():"*")+"="+getPeptide()+"^"+
         MsStats.asDeviation(mDeltaM, mCalcMH, 999) + "->" + Tools.d2s(getGapScore(), 2);
  }

  @Override
  public int compareTo(Ms2Hit o)
  {
    int c = Long.compare(mProteinKey, o.getProteinKey());

    if (c==0) c = Integer.compare(mLeft,  o.getLeft());
    if (c==0) c = Integer.compare(mRight, o.getRight());
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
        if (tol.withinTolerance(getDelta()+getCalcMH(), i*1.00335d+getCalcMH()))
        {
          mIsotopeError=i; mDeltaM-=(i*1.00335d);
          break;
        }

    return this;
  }
  @Override
  public void dispose()
  {
    Tools.dispose(mY, mB);
    mProteinKey=null;
    mSequence=mPrev=mNext=null;
    Tools.dispose(mMods);
  }
}
