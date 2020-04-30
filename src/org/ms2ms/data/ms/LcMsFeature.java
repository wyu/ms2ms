package org.ms2ms.data.ms;

import org.ms2ms.data.Point;

/** A representation of LcMs feature from MaxQuant or similar program.
 *
 */
public class LcMsFeature extends LcMsPoint
{
  private double mArea;
  private double mAbundance;
  private double mApex;

  public LcMsFeature() {};
  public LcMsFeature(double mz, double rt)
  {
  }
  public LcMsFeature(Point s)
  {
    if (s!=null) { setX(s.getX()); setY(s.getY()); }
  }

  public double getAbundance() { return mAbundance; }

  public LcMsPoint setAbundance(double s) { mAbundance=s; return this; }

  public double getArea()      { return mArea; }

  public LcMsPoint setArea(double s) { mArea=s; return this; }

  public double getApex()      { return mApex; }

  public LcMsPoint setApex(double s) { mApex=s; return this; }
}
