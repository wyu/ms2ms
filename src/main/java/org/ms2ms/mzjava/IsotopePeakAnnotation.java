package org.ms2ms.mzjava;

import org.expasy.mzjava.core.ms.peaklist.PeakAnnotation;

/**
 * Created by yuw on 8/9/16.
 */
public class IsotopePeakAnnotation implements PeakAnnotation
{
  private int mIsotopeOrder=0, mCharge=0;
  private double mIntensity=0d;

  public IsotopePeakAnnotation() { super(); }
  public IsotopePeakAnnotation(IsotopePeakAnnotation s)
  {
    super();
    mCharge=s.mCharge; mIsotopeOrder=s.mIsotopeOrder; mIntensity=s.mIntensity;
  }
//  public IsotopePeakAnnotation(int z, int s)
//  {
//    super();
//    mCharge=z; mIsotopeOrder=s;
//  }
  public IsotopePeakAnnotation(int z, int s, double ai)
  {
    super();
    mCharge=z; mIsotopeOrder=s; mIntensity=ai;
  }

  public int getIsotopeOrder() { return mIsotopeOrder; }

  @Override
  public int getCharge() { return mCharge; }

  public double getIntensity() { return mIntensity; }
  public IsotopePeakAnnotation setIntensity(double s) { mIntensity=s; return this; }

  @Override
  public PeakAnnotation copy() { return new IsotopePeakAnnotation(this); }

  @Override
  public boolean equals(Object o) {

    if (this == o) return true;
    if (!(o instanceof IsotopePeakAnnotation)) return false;

    IsotopePeakAnnotation that = (IsotopePeakAnnotation) o;

    return that.mCharge==mCharge && that.mIsotopeOrder==mIsotopeOrder;
  }
  @Override
  public int hashCode()
  {
    return Integer.valueOf(mIsotopeOrder+mCharge+mIsotopeOrder*mCharge).hashCode();
  }
  @Override
  public String toString() { return "z"+mCharge+" @"+mIsotopeOrder+"#"+mIntensity; }
}
