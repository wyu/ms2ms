package org.ms2ms.splib;

import com.google.common.collect.Range;
import com.google.common.collect.TreeMultimap;
import org.apache.commons.math3.analysis.BivariateFunction;
import org.apache.commons.math3.analysis.interpolation.PiecewiseBicubicSplineInterpolator;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.utils.Tools;

import java.util.Collections;

public class Sage2D extends AbstractSage
{
  private Integer       mBinX=25;
  private Range<Double> mBound, mBoundX;

  private MultiTreeTable<Double, Double, Integer> mPositives, mNegatives;
  private BivariateFunction mTransition;

  public Sage2D() { super(); }
  public Sage2D(String s, int bins, int binsX) { mTitle=s; mBins=bins; mBinX=binsX; }

  public void addPositive(double data, double x)
  {
    if (mPositives==null) mPositives = MultiTreeTable.create();
    mPositives.put(data,x,1);
    mPosCounts++;
  }
  public void addNegative(double data, double x)
  {
    if (mNegatives==null) mNegatives = MultiTreeTable.create();
    mNegatives.put(data,x,1);
    mNegCounts++;
  }

  private void updateBounds()
  {
    if (mBound==null) mBound = Range.closed(
        Math.min(Collections.min(mPositives.keySet()), Collections.min(mNegatives.keySet())),
        Math.max(Collections.max(mPositives.keySet()), Collections.max(mNegatives.keySet())));
    if (mBoundX==null) mBoundX = Range.closed(
        Math.min(Collections.min(mPositives.columnSet()), Collections.min(mNegatives.columnSet())),
        Math.max(Collections.max(mPositives.columnSet()), Collections.max(mNegatives.columnSet())));
  }
  protected void toTransition(double p_pos, double p_neg, boolean keep_zero)
  {
    if (mPositives == null || mNegatives == null)
      throw new RuntimeException("The positive and/or negative population not specified");

    mPosPrior=p_pos; mNegPrior=p_neg;

    updateBounds();

    double step = ( mBound.upperEndpoint()- mBound.lowerEndpoint()) / (mBins*2d),
          xstep = (mBoundX.upperEndpoint()-mBoundX.lowerEndpoint()) / (mBinX*2d),
         lowest = Double.MAX_VALUE, base_pos = mPosCounts/(mBins*2d), base_neg = mNegCounts/(mBins*2d);

    double[] pts = new double[mBins*2-1], ptsX = new double[mBinX*2-1];
    double[][] pp = new double[mBins*2-1][mBinX*2-1];

    for (  int i=0; i<mBins*2-1; i++)
      for (int j=0; j<mBinX*2-1; j++)
      {
        pts[i] = mBound.lowerEndpoint()+i*step; ptsX[j] = mBoundX.lowerEndpoint()+j*xstep;
        double pd_pos = mPositives.subset(pts[i],pts[i]+step, ptsX[j],ptsX[j]+xstep).size() / base_pos,
               pd_neg = mNegatives.subset(pts[i],pts[i]+step, ptsX[j],ptsX[j]+xstep).size() / base_neg;

      pd_pos = pd_pos < 0 ? 0 : pd_pos;
      pd_neg = pd_neg < 0 ? 0 : pd_neg;

      double p = (pd_pos == 0 && pd_neg == 0) ? 0d : pd_pos * p_pos / (pd_pos * p_pos + pd_neg * p_neg);

      pp[i][j] = !Double.isNaN(p)?p:0d;
      if (p!=0 && p<lowest) lowest=p;
    }
    if (!keep_zero)
      for (  int i=0; i<mBins*2-1; i++)
        for (int j=0; j<mBinX*2-1; j++)
          if (pp[i][j]==0) pp[i][j] = lowest*0.1;

    mTransition = new PiecewiseBicubicSplineInterpolator().interpolate(pts, ptsX, pp);
  }
}
