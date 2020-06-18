package org.ms2ms.data.ms;

import org.ms2ms.data.Point;

/** A representation of LcMs feature from MaxQuant or similar program.
 *
 */
public class LcMsFeature extends LcMsPoint
{
  private double mArea;
  private double mAbundance;
  private double mApex, mDeriv1st;
  private double mMzStdev=Double.NaN, mMzStdevEx=Double.NaN;
  private double mInitialCentroidRt;

  public LcMsFeature() {};
  public LcMsFeature(double mz, double rt)
  {
  }
  public LcMsFeature(double rt, double ai, double d1)
  {
    setX(rt); setY(ai); mDeriv1st=d1;
  }
  public LcMsFeature(Point s)
  {
    if (s!=null) { setX(s.getX()); setY(s.getY()); }
  }

  public double getAbundance() { return mAbundance; }
  public double getArea()      { return mArea; }
  public double getApex()      { return mApex; }
  public double getMzStdev()   { return mMzStdev; }
  public double getMzStdevEx() { return mMzStdevEx; }
  public double getInitialCentroidRt()      { return mInitialCentroidRt; }

  public LcMsPoint setAbundance(double s) { mAbundance=s; return this; }
  public LcMsPoint setArea(double s) { mArea=s; return this; }
  public LcMsPoint setApex(double s) { mApex=s; return this; }
//  public LcMsPoint setMz(double s) { mMz=s; return this; }
  public LcMsPoint setMzStdev(double s) { mMzStdev=s; return this; }
  public LcMsPoint setMzStdevEx(double s) { mMzStdevEx=s; return this; }

  public LcMsFeature setInitialCentroidRt(double s) { mInitialCentroidRt=s; return this; }

}
