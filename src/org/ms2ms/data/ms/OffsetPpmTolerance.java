package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.PpmTolerance;

/**
 * Created by yuw on 8/6/16.
 */
public class OffsetPpmTolerance extends PpmTolerance
{
  private double mScale=1d, mOffset, mTol;

  public OffsetPpmTolerance() { super(0d); }
  public OffsetPpmTolerance(double tol, double offset)
  {
    super(tol); mOffset=offset; mTol=tol;
  }

  public OffsetPpmTolerance scale(double s) { mScale=s; return this; }

  @Override
  public Location check(double expected, double actual)
  {
    if      (actual < getMin(expected)) return Location.SMALLER;
    else if (actual > getMax(expected)) return Location.LARGER;
    else                                return Location.WITHIN;
  }

  @Override
  public boolean withinTolerance(double expected, double actual) { return actual >= getMin(expected) && actual <= getMax(expected); }

  @Override
  public double getMin(double mz) { return  mz-calcError(mz)+calcOffset(mz); }
  @Override
  public double getMax(double mz) { return  mz+calcError(mz)+calcOffset(mz); }

  private double calcError( double expectedMass) { return expectedMass * (mTol*mScale/1000000d); }
  private double calcOffset(double expectedMass) { return expectedMass * (mOffset/1000000d); }
}
