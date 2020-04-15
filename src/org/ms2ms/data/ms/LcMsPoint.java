package org.ms2ms.data.ms;

import org.ms2ms.data.Point;

public class LcMsPoint extends Point
{
  private double mMz;
  private int mScan;

  public LcMsPoint() { super(); }
  public LcMsPoint(double rt, double ai) { super(rt, ai); }
  public LcMsPoint(double rt, double ai, double mz, int scan)
  {
    super(rt, ai);
    mMz=mz; mScan=scan;
  }

  public double getMz()        { return mMz; }
  public double getRT()        { return getX(); }
  public double getIntensity() { return getY(); }
  public int    getScan()      { return mScan; }

  public LcMsPoint setMz(       double s) { mMz  =s;  return this; }
  public LcMsPoint setScan(        int s) { mScan=s;  return this; }
  public LcMsPoint setRT(       double s) { setX( s); return this; }
  public LcMsPoint setIntensity(double s) { setY( s); return this; }
}
