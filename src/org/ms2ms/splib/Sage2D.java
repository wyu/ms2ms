package org.ms2ms.splib;

import com.google.common.collect.Range;
import com.google.common.collect.TreeBasedTable;
import com.google.common.collect.TreeMultimap;
import org.apache.commons.math3.analysis.BivariateFunction;
import org.apache.commons.math3.analysis.interpolation.PiecewiseBicubicSplineInterpolator;
import org.ms2ms.algo.MsStats;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.utils.Tools;

import java.util.Collections;

public class Sage2D extends AbstractSage
{
  private Integer       mBinX=25, mSpan=4;
  private Range<Double> mBound, mBoundX;

  private TreeBasedTable<Double, Double, Integer> mPositives, mNegatives;
  private BivariateFunction mTransition;

  public Sage2D() { super(); }
  public Sage2D(String s, int bins, int binsX) { mTitle=s; mBins=bins; mBinX=binsX; }

  public void addPositive(double data, double x)
  {
    mPositives = add(mPositives, data,x);
    mPosCounts++;
  }
  public void add(double data, double x, boolean positive) { if (positive) addPositive(data,x); else addNegative(data,x); }
  public void addNegative(double data, double x)
  {
    mNegatives = add(mNegatives, data,x);
    mNegCounts++;
  }
  private TreeBasedTable<Double, Double, Integer> add(TreeBasedTable<Double, Double, Integer> pts, double data, double x)
  {
    if (pts==null) pts = TreeBasedTable.create();
    pts.put(data,x,pts.get(data,x)!=null?pts.get(data,x)+1:1);
    return pts;
  }

  private void updateBounds()
  {
    if (mBound==null) mBound = Range.closed(
        Math.min(Collections.min(mPositives.rowKeySet()), Collections.min(mNegatives.rowKeySet())),
        Math.max(Collections.max(mPositives.rowKeySet()), Collections.max(mNegatives.rowKeySet())));
    if (mBoundX==null) mBoundX = Range.closed(
        Math.min(Collections.min(mPositives.columnKeySet()), Collections.min(mNegatives.columnKeySet())),
        Math.max(Collections.max(mPositives.columnKeySet()), Collections.max(mNegatives.columnKeySet())));
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
        Range<Double> r1 = Tools.bound(mBound,  pts[i] -mSpan* step, pts[i] +mSpan* step),
                      r2 = Tools.bound(mBoundX, ptsX[j]-mSpan*xstep, ptsX[j]+mSpan*xstep);
        double pd_pos = Tools.sliceCounts(mPositives, r1, r2) / base_pos,
               pd_neg = Tools.sliceCounts(mNegatives, r1, r2) / base_neg;

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

  public double lookup(double data, double x) { return lookup(data,x,mPosPrior, mNegPrior); }
  public double lookup(double data, double x, double p_pos, double p_neg)
  {
    if (mPositives == null || mNegatives == null)
      throw new RuntimeException("The positive and/or negative population not specified");

    if (mTransition==null) toTransition(p_pos, p_neg, false);

    return mTransition.value(data,x); // ignore zero?
  }

}
