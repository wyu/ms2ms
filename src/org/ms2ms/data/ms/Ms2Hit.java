package org.ms2ms.data.ms;

import org.ms2ms.algo.MsStats;
import org.ms2ms.utils.Tools;

/**
 * Created by yuw on 8/7/16.
 */
public class Ms2Hit implements Comparable<Ms2Hit>
{
  private FpmEntry mY, mB;
  private Long mProteinKey;
  private int mLeft, mRight;
  private double mCalcMH, mDeltaM, mProb;
  private String mPeptide;

  public Ms2Hit()  { super(); }
  public Ms2Hit(Long protein, FpmEntry y, FpmEntry b, int left, int right)
  {
    super();
    mProteinKey=protein; mY=y; mB=b; mLeft=left; mRight=right;
  }

  public Long     getProteinKey() { return mProteinKey; }
  public int      getLeft()       { return mLeft; }
  public int      getRight()      { return mRight; }
  public double   getDelta()      { return mDeltaM; }
  public double   getCalcMH()     { return mCalcMH; }
  public double   getProb()       { return mProb; }
  public FpmEntry getY()          { return mY; }
  public FpmEntry getB()          { return mB; }
  public int      getMotifs()     { return(mY!=null?mY.getMotifs():0)+(mB!=null?mB.getMotifs():0); }
  public String   getPeptide()    { return mPeptide; }

  public Ms2Hit   setLeft( int s) { mLeft =s; return this; }
  public Ms2Hit   setRight(int s) { mRight=s; return this; }
  public Ms2Hit   setPeptide(String s) { mPeptide=s; return this; }
  public Ms2Hit   setProb(double s) { mProb=s; return this; }
  public Ms2Hit   setMH(double calc, double delta)
  {
    mCalcMH=calc; mDeltaM=delta;
    return this;
  }
  public Ms2Hit setProb()
  {
    mProb=(getY()!=null?getY().getProb():0d)+(getB()!=null?getB().getProb():0d);
    return this;
  }
  @Override
  public String toString()
  {
    return getProteinKey()+":"+getLeft()+"-"+getRight()+",m/z"+ Tools.d2s(mCalcMH, 5)+"$"+
        (mB!=null?mB.getTrack().size():"*")+"->"+(mY!=null?mY.getTrack().size():"*")+"="+getPeptide()+"^"+ MsStats.asDeviation(mDeltaM, mCalcMH);
  }

  @Override
  public int compareTo(Ms2Hit o)
  {
    int c = Long.compare(mProteinKey, o.getProteinKey());

    if (c==0) c = Integer.compare(mLeft,  o.getLeft());
    if (c==0) c = Integer.compare(mRight, o.getRight());

    return c;
  }
}
