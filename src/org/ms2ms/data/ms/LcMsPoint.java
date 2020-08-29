package org.ms2ms.data.ms;

import com.google.common.collect.Range;
import org.ms2ms.data.Point;
import org.ms2ms.math.Points;
import org.ms2ms.utils.Tools;

import java.io.IOException;
import java.io.Writer;
import java.util.Collection;
import java.util.List;

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
  public static LcMsPoint interpolate(List<LcMsPoint> ps, Double x, boolean ignore_zero)
  {
    Range<LcMsPoint> range = Points.boundry(new LcMsPoint(x, 0d), ps, ignore_zero);
    if (range!=null)
      return interpolate(range.lowerEndpoint(), range.upperEndpoint(), x);

    return null;
  }
  public static LcMsPoint interpolate(LcMsPoint p1, LcMsPoint p2, Double x)
  {
    LcMsPoint xy = new LcMsPoint(x, 0d);

    if     ((p1 != null && p2 == null) || (p1 == null && p2 != null)) xy = null; // undefined situation, WYU 081209
    else if (p1 != null && p2 != null && p2.getX() - p1.getX() != 0)
    {
      Double k0 = (p2.getX() - p1.getX()),
              k = (p2.getY()        - p1.getY())        / k0,
            fil = (p2.getFillTime() - p1.getFillTime()) / k0,
            ppm = (p2.getPPM()      - p1.getPPM())      / k0,
             mz = (p2.getMz()       - p1.getMz())       / k0,
           scan = (p2.getScan()     - p1.getScan())     / k0;

      xy.setY(       p1.getY()        + (x - p1.getX()) * k);
      xy.setFillTime(p1.getFillTime() + (x - p1.getX()) * fil);
      xy.setMz(      p1.getMz()       + (x - p1.getX()) * mz);
      xy.setPPM((float          )(p1.getPPM()  + (x - p1.getX()) * ppm));
      xy.setScan((int )Math.round(p1.getScan() + (x - p1.getX()) * scan));
    }
    else if (p1 != null && p2 != null && p2.getX() == p1.getX())
    {
      xy.setY(       0.5d*(p1.getY() + p2.getY()));
      xy.setFillTime(0.5d*(p1.getFillTime() + p2.getFillTime()));
      xy.setMz(      0.5d*(p1.getMz() + p2.getMz()));
      xy.setPPM((float          )(0.5f*(p1.getPPM()  + p2.getPPM())));
      xy.setScan((int )Math.round(0.5f*(p1.getScan() + p2.getScan())));
    }
    return xy;
  }

  @Override
  public String toString()
  {
    String out = "rt "+Tools.d2s(getRT(), 2) + ", m/z " + Tools.d2s(getMz(), 4) + "," +
        Tools.d2s(getIntensity(), 1) + "," + Tools.d2s(getFillTime(), 1) + " msec";
    return out;
  }
}
