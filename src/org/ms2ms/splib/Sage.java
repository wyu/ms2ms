package org.ms2ms.splib;

import com.google.common.collect.Range;
import org.expasy.mzjava.stats.Histogram;
import org.ms2ms.algo.MsStats;
import org.ms2ms.utils.Tools;

/** Sage: Empirical Bayes?
 *
 * Date: 8/9/14
 * Time: 8:59 AM
 */
public class Sage
{
    private String        mTitle;
    private int           mBins = 50;
    private Range<Double> mBound;
    private Histogram     mPositives, mNegatives;
    private double[]      mTransitXs, mTransitYs;

    public Sage() { super(); }
    public Sage(String title, Range<Double> bound, int bins)
    {
      super();
      mTitle=title; mBins=bins; mBound=bound;
      mPositives= MsStats.newHistogram("Positives", mBins, mBound);
      mNegatives= MsStats.newHistogram("Negatives", mBins, mBound);
    }

    public String      getTitle()      { return mTitle; }
    public Histogram   getPositives()  { return mPositives; }
    public Histogram   getNegatives()  { return mNegatives; }
    public double[]    getTransitXs()  { return mTransitXs; }
    public double[]    getTransitYs()  { return mTransitYs; }

    public void addPositive(double data)
    {
      if (mPositives == null) mPositives= MsStats.newHistogram("Positives", mBins, mBound);
      mPositives.addData(data);
    }
    public void addNegative(double data)
    {
      if (mNegatives == null) mNegatives= MsStats.newHistogram("Negatives", mBins, mBound);
      mNegatives.addData(data);
    }

    public double lookup(double x, double p_pos, double p_neg)
    {
      if (mPositives == null || mNegatives == null)
        throw new RuntimeException("The positive and/or negative population not specified");

      if (mTransitYs==null || mTransitXs==null) toTransition(p_pos, p_neg);

      if      (x>=Tools.back( mTransitXs)) return Tools.back( mTransitYs);
      else if (x<=Tools.front(mTransitXs)) return Tools.front(mTransitYs);

      return MsStats.interpolate(mTransitXs, mTransitYs, 0.5d, x)[0]; // ignore zero?
    }
    private void toTransition(double p_pos, double p_neg)
    {
      if (mPositives == null || mNegatives == null)
        throw new RuntimeException("The positive and/or negative population not specified");

      mTransitXs=new double[mBins*2]; mTransitYs=new double[mBins*2];

      // normalize the frequency to the total area of 1
      mPositives.normalize(new Histogram.Normalization(Histogram.Normalization.NormType.BINS_CUMUL, 1d));
      mNegatives.normalize(new Histogram.Normalization(Histogram.Normalization.NormType.BINS_CUMUL, 1d));

      double xstep = (mBound.upperEndpoint()-mBound.lowerEndpoint()) / (mBins*2d); int i=0;
      for (double x=mBound.lowerEndpoint(); x <= mBound.upperEndpoint(); x += xstep)
      {
        double pd_pos = mPositives.getAbsoluteBinFreq(mPositives.getBinIndex(x)),
               pd_neg = mNegatives.getAbsoluteBinFreq(mNegatives.getBinIndex(x));

        pd_pos = pd_pos < 0 ? 0 : pd_pos;
        pd_neg = pd_neg < 0 ? 0 : pd_neg;

        double p = (pd_pos == 0 && pd_neg == 0) ? 0d : pd_pos * p_pos / (pd_pos * p_pos + pd_neg * p_neg);

        if (!Double.isNaN(p))
        {
          mTransitXs[i]=x; mTransitYs[i]=p;
        }
        i++;
      }
    }
}
