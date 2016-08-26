package org.ms2ms.data.ms;

import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
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
  private FpmEntry mY, mB;
  private Long mProteinKey;
  private int mLeft, mRight, mRank;
  private double mCalcMH, mDeltaM, mKaiScore, mDeltaScore, mEval, mScoreOffset=0;
  private String mSequence, mPrev, mNext;
  private TreeMap<Integer, Double> mMods;

  public Ms2Hit()
  {
    super();
    mY = new FpmEntry(); mB = new FpmEntry();
  }
  public Ms2Hit(Long protein, FpmEntry y, FpmEntry b, int left, int right)
  {
    super();
    mProteinKey=protein; mY=y; mB=b; mLeft=left; mRight=right;
  }

  public boolean  isDecoy()       { return mProteinKey!=null && mProteinKey<0; }
  public Long     getProteinKey() { return mProteinKey; }
  public int      getLeft()       { return mLeft; }
  public int      getRight()      { return mRight; }
  public int      getRank()       { return mRank; }
  public double   getDelta()      { return mDeltaM; }
  public double   getCalcMH()     { return mCalcMH; }
  public double   getEVal()       { return mEval; }
  public double   getKaiScore()   { return mKaiScore; }
  public double   getDeltaScore() { return mDeltaScore; }
  public FpmEntry getY()          { return mY; }
  public FpmEntry getB()          { return mB; }
//  public int      getMotifs()     { return(mY!=null?mY.getMotifs(  ):0)+(mB!=null?mB.getMotifs(  ):0); }
  public double   getGapScore()   { return(mY!=null?mY.getGapScore():0)+(mB!=null?mB.getGapScore():0)+mScoreOffset; }
  public String   getPeptide()    { return (mPrev+"."+mSequence+"."+mNext); }
  public String   getSequence()   { return mSequence; }
  public double   getModMass()    { return mMods!=null? Stats.sum(mMods.values()):0d; }
//  public TreeMap<Integer, Double> getMods() { return mMods; }

  public Map<Integer, Double> getMod0()
  {
    if (mMods==null) return null;

    SortedMap<Integer, Double> out = new TreeMap<>();
    for (Map.Entry<Integer, Double> E : mMods.entrySet())
      out.put(E.getKey()-getLeft(), E.getValue());

    return out;
  }

  public Ms2Hit   setLeft(          int s) { mLeft =s; return this; }
  public Ms2Hit   setRight(         int s) { mRight=s; return this; }
  public Ms2Hit   setRank(          int s) { mRank=s; return this; }
  //public Ms2Hit   setPeptide(String s) { mPeptide=s; return this; }
  public Ms2Hit   setEVal(       double s) { mEval =s; return this; }
  public Ms2Hit   setKaiScore(   double s) { mKaiScore=s; return this; }
  public Ms2Hit   setDeltaScore( double s) { mDeltaScore=s; return this; }
  public Ms2Hit   setScoreOffset(double s) { mScoreOffset=s; return this; }
  public Ms2Hit   setMH(double calc, double delta)
  {
    mCalcMH=calc; mDeltaM=delta;
    return this;
  }
  public Ms2Hit setPeptide(String sequence) { return sequence!=null?setPeptide(sequence.toCharArray(), 0, sequence.length()-1):this; }
  public Ms2Hit setPeptide(char[] sequence) { return setPeptide(sequence, getLeft(), getRight()); }
  public Ms2Hit setPeptide(char[] sequence, int left, int right)
  {
    setLeft(left); setRight(right);
    mPrev    =(getLeft()>0?sequence[getLeft()-1]:'-')+"";
    mSequence=Strs.toString(sequence, getLeft(), getRight()+1);
    mNext    =(getRight()<sequence.length-1?sequence[getRight()+1]:'-')+"";

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
      throw new RuntimeException("Conflicting mod from N/C terminus!");
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
         MsStats.asDeviation(mDeltaM, mCalcMH) + "->" + Tools.d2s(getGapScore(), 2);
  }

  @Override
  public int compareTo(Ms2Hit o)
  {
    int c = Long.compare(mProteinKey, o.getProteinKey());

    if (c==0) c = Integer.compare(mLeft,  o.getLeft());
    if (c==0) c = Integer.compare(mRight, o.getRight());

    return c;
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
