package org.ms2ms.splib;

abstract public class AbstractSage
{
  protected String        mTitle;
  protected int           mBins = 50;
  protected int mPosCounts=0;
  protected int mNegCounts=0;
  protected double        mPosPrior=0.5d;
  protected double mNegPrior=0.5d;

  public String      getTitle()      { return mTitle; }

  public void setPriors(double p_pos, double p_neg) { mPosPrior=p_pos; mNegPrior=p_neg; }
}
