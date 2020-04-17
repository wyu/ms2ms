package org.ms2ms.data.ms;

import org.ms2ms.data.Point;
import org.ms2ms.utils.Tools;

public class LcMsPoint extends Point
{
  private double  mMz=0, mApex, mArea, mAbundance, mFillTime;
  private int     mScan=0;
  private boolean mImputed=false;

  public LcMsPoint() { super(); }
  public LcMsPoint(double rt, double ai) { super(rt, ai); }
  public LcMsPoint(double rt, double ai, double mz, int scan)
  {
    super(rt, ai);
    mMz=mz; mScan=scan;
  }
  public LcMsPoint(Point s)
  {
    if (s!=null) { setX(s.getX()); setY(s.getY()); }
  }

  public double getMz()        { return mMz; }
  public double getRT()        { return getX(); }
  public double getApex()      { return mApex; }
  public double getArea()      { return mArea; }
  public double getAbundance() { return mAbundance; }
  public double getIntensity() { return getY(); }
  public double getFillTime()  { return mFillTime; }
  public int    getScan()      { return mScan; }
  public boolean isImputed()   { return mImputed; }

  public LcMsPoint setMz(       double s) { mMz  =s;  return this; }
  public LcMsPoint setScan(        int s) { mScan=s;  return this; }
  public LcMsPoint setRT(       double s) { setX( s); return this; }
  public LcMsPoint setApex(     double s) { mApex=s; return this; }
  public LcMsPoint setArea(     double s) { mArea=s; return this; }
  public LcMsPoint setAbundance(double s) { mAbundance=s; return this; }
  public LcMsPoint setIntensity(double s) { setY( s); return this; }
  public LcMsPoint setFillTime( double s) { mFillTime=s; return this; }
  public LcMsPoint isImputed(  boolean s) { mImputed=s; return this; }

  @Override
  public String toString()
  {
    String out = "rt:"+Tools.d2s(getRT(), 2) + ", m/z" + Tools.d2s(getMz(), 4) + "->" +
        Tools.d2s(getIntensity(), 1) + ", " + Tools.d2s(getFillTime(), 1) + "msec";
    return out;
  }
}
