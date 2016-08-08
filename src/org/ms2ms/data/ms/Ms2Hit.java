package org.ms2ms.data.ms;

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

  public Long   getProteinKey() { return mProteinKey; }
  public int    getLeft()       { return mLeft; }
  public int    getRight()      { return mRight; }
  public double getDelta()      { return mDeltaM; }
  public double getCalcMH()     { return mCalcMH; }

  public Ms2Hit setLeft( int s) { mLeft =s; return this; }
  public Ms2Hit setRight(int s) { mRight=s; return this; }
  public Ms2Hit setMH(double calc, double delta)
  {
    mCalcMH=calc; mDeltaM=delta;
    return this;
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
