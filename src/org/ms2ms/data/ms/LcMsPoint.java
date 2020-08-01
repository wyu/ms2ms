package org.ms2ms.data.ms;

import org.ms2ms.data.Point;
import org.ms2ms.utils.Tools;

import java.io.IOException;
import java.io.Writer;
import java.util.Collection;

public class LcMsPoint extends Point
{
  private float  mMz=0, mPPM=Float.NaN;
  private float  mFillTime;
  private int     mScan=0;
  private boolean mImputed=false;

  public LcMsPoint() { super(); }
  public LcMsPoint(Point s) { super(s.getX(), s.getY()); }
  public LcMsPoint(double rt, double ai)
  {
    super(rt, ai);
  }
  public LcMsPoint(double rt, double ai, float mz, int scan, float ppm)
  {
    super(rt, ai);
    mMz=mz; mScan=scan; mPPM=ppm;
  }

  public double getMz()        { return mMz; }
  public double getRT()        { return getX(); }
  public double getPPM()       { return mPPM; }

  public double getIntensity() { return getY(); }
  public double getFillTime()  { return mFillTime; }
  public int    getScan()      { return mScan; }
  public boolean isImputed()   { return mImputed; }

  public LcMsPoint setPPM(      float s) { mPPM  =s;  return this; }
  public LcMsPoint setMz(       Double s) { mMz  =(s!=null?(float )s.doubleValue():0f);  return this; }
  public LcMsPoint setMz(Double s, float f) { setMz(s); mPPM = (s!=null && f!=0)?(float )(1E6*(s-f)/f):Float.NaN;  return this; }

  public LcMsPoint setScan(        int s) { mScan=s;  return this; }
  public LcMsPoint setRT(       double s) { setX( s); return this; }

  public LcMsPoint setFillTime(  Double s) { mFillTime=(s!=null?(float )s.doubleValue():0f); return this; }
  public LcMsPoint isImputed(  boolean s) { mImputed=s; return this; }
  public LcMsPoint setIntensity(double s)
  {
    if (Double.isNaN(s) || Double.isInfinite(s) || Double.isNaN(getY()) || Double.isInfinite(getY()))
      System.out.println();
    setY( s); return this;
  }
  public static <T extends LcMsPoint> Double sumY(Collection<T> data)
  {
    if (data==null) return null;

    double sum = 0d;
    for (T xy : data) sum += xy.getY();

    return sum;
  }
  public void printXIC(Writer w) throws IOException
  {
    w.write(Tools.d2s(getScan(),2)     +"\t");
    w.write(Tools.d2s(getRT(),3)       +"\t");
    w.write(Tools.d2s(getIntensity(),2)+"\t");
    w.write(Tools.d2s(getFillTime(),2)+"\t");
    w.write(isImputed()                   +"\t");
    w.write(Tools.d2s(getPPM(),2)      +"\n");
  }

  @Override
  public String toString()
  {
    String out = "rt "+Tools.d2s(getRT(), 2) + ", m/z " + Tools.d2s(getMz(), 4) + "," +
        Tools.d2s(getIntensity(), 1) + "," + Tools.d2s(getFillTime(), 1) + " msec";
    return out;
  }
}
