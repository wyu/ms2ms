package org.ms2ms.data.ms;

import org.ms2ms.utils.Tools;

/**
 * Created by yuw on 8/7/16.
 */
public class Ms2Hit implements Comparable<Ms2Hit>
{
  private FpmEntry mY, mB;
  private Long mProteinKey;
  private int mLeft, mRight;
  private double mCalcMH, mDeltaM;

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
  public FpmEntry getY()          { return mY; }
  public FpmEntry getB()          { return mB; }

  public Ms2Hit setLeft( int s) { mLeft =s; return this; }
  public Ms2Hit setRight(int s) { mRight=s; return this; }
  public Ms2Hit setMH(double calc, double delta)
  {
    mCalcMH=calc; mDeltaM=delta;
    return this;
  }
  @Override
  public String toString()
  {
    return getProteinKey()+":"+getLeft()+"-"+getRight()+",m/z"+ Tools.d2s(mCalcMH, 5)+"$"+
        (mB!=null?mB.getTrack().size():"*")+"->"+(mY!=null?mY.getTrack().size():"*");
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
