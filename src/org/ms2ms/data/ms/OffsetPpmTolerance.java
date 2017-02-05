package org.ms2ms.data.ms;

import com.google.common.collect.Range;
import org.expasy.mzjava.core.ms.PpmTolerance;
import org.ms2ms.utils.Tools;

/**
 * Created by yuw on 8/6/16.
 */
public class OffsetPpmTolerance extends PpmTolerance
{
  private double mScale=1d, mOffset=0d, mTol=0d, mTransitionMass=0, mOffsetSlope=0, mTolSlope=0, mZval=0;
  private boolean mIsIncremental=false;

  public OffsetPpmTolerance() { super(0d); }
  public OffsetPpmTolerance(double tol) { super(tol); mTol=tol; }
  public OffsetPpmTolerance(double tol, double offset)
  {
    super(tol); mOffset=offset; mTol=tol;
  }

  public OffsetPpmTolerance isIncremental(boolean s) { mIsIncremental=s; return this; }
  public boolean            isIncremental()          { return mIsIncremental; }

  public OffsetPpmTolerance setTransitionMass(double s) { mTransitionMass=s; return this; }
  public OffsetPpmTolerance setOffsetParams(double transition, double intercept, double slope)
  {
    mTransitionMass=transition; mOffsetSlope=slope; mOffset=intercept; return this;
  }
  public OffsetPpmTolerance setTolParams(double intercept, double slope, double zval)
  {
    mTolSlope=slope; mTol=intercept; mZval=zval; return this;
  }
  public double getOffset(double m)
  {
    // flip the sign to indicate correction
    if (!isIncremental() && mOffsetSlope!=0) return m<mTransitionMass?(-mOffset-mOffsetSlope*m):(getOffset(mTransitionMass-1d));
    return mOffset;
  }
  // 2-sigma to cover the 95% interval
  public double getPpmTol(double m)
  {
    if (isIncremental())
    {
      return m<mTransitionMass?(mTol+mTolSlope*m)*mZval*mScale:(getPpmTol(mTransitionMass-1d));
    }
    return mTolSlope!=0?Math.exp(mTol+mTolSlope*m)*mZval*mScale : (mTol*mScale);
  }
//  public double getPpmTol() { return mTol; }

  public OffsetPpmTolerance scale( double s) { mScale =s; return this; }
  public OffsetPpmTolerance offset(double s) { mOffset=s; return this; }

  @Override
  public Location check(double expected, double actual)
  {
    if      (actual < getMin(expected)) return Location.SMALLER;
    else if (actual > getMax(expected)) return Location.LARGER;
    else                                return Location.WITHIN;
  }

//  @Deprecated
  public Range<Double> getBoundary(double mz) { return Range.closed(getMin(mz), getMax(mz)); }
  public double[] toActualBoundary(  double mz) { return new double[] {mz-calcError(mz)-calcOffset(mz), mz+calcError(mz)-calcOffset(mz)}; }
  public float[]  toExpectedBoundary(float mz)  { return new float[] {(float )(mz-calcError(mz)+calcOffset(mz)), (float )(mz+calcError(mz)+calcOffset(mz))}; }

  @Override
  public boolean withinTolerance(double expected, double actual) { return actual >= getMin(expected) && actual <= getMax(expected); }

  @Override
  public double getMin(double mz) { return  mz-calcError(mz)+calcOffset(mz); }
  @Override
  public double getMax(double mz) { return  mz+calcError(mz)+calcOffset(mz); }

  public OffsetPpmTolerance clone()
  {
    return new OffsetPpmTolerance().setOffsetParams(mTransitionMass, mOffset, mOffsetSlope).setTolParams(mTol, mTolSlope, mZval).isIncremental(mIsIncremental);
  }
  private double calcError( double expectedMass) { return expectedMass * (getPpmTol(expectedMass)/1000000d); }
  private double calcOffset(double expectedMass) { return expectedMass * (getOffset(expectedMass)/1000000d); }

  @Override
  public String toString()
  {
    String scale = (mScale==1?("*"+Tools.d2s(mScale,2)):"");
    String ppm = (mTolSlope   !=0 ? "exp("+Tools.d2s(mTol,2)+"+"+Tools.d2s(mTolSlope,7)+"*m)*"+Tools.d2s(mZval,2)+scale : (Tools.d2s(mTol,2)+scale));
    String oft = (mOffsetSlope!=0 ? (Tools.d2s(-mOffset,2)+"-"+Tools.d2s(mOffsetSlope,7)+"*m or " + Tools.d2s(getOffset(mTransitionMass-1D),2) + " above m/z") : "") + Tools.d2s(mTransitionMass, 2);

    return "ppm = "+ppm+ ", offset = "+oft+ (mIsIncremental?", Incremental":"");
  }
}
