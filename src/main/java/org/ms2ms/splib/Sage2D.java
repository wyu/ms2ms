package org.ms2ms.splib;

import com.google.common.collect.Range;
import com.google.common.collect.TreeBasedTable;
import org.apache.commons.math3.analysis.BivariateFunction;
import org.apache.commons.math3.analysis.interpolation.PiecewiseBicubicSplineInterpolator;
import org.apache.commons.math3.exception.OutOfRangeException;
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
    if (!Double.isInfinite(x) && !Double.isNaN(x) && !Double.isInfinite(data) && !Double.isNaN(data))
    {
      if (pts==null) pts = TreeBasedTable.create();
      pts.put(data,x,pts.get(data,x)!=null?pts.get(data,x)+1:1);
    }
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

    double[] pts = new double[mBins*2-1], ptsX = new double[mBinX*2-1], multis = new double[]{1d,2d,4d,8d,16d};
    double[][] pp = new double[mBins*2-1][mBinX*2-1];

    int minN=16;
    for (  int i=0; i<mBins*2-1; i++)
      for (int j=0; j<mBinX*2-1; j++)
      {
        pts[i] = mBound.lowerEndpoint()+i*step; ptsX[j] = mBoundX.lowerEndpoint()+j*xstep;
        double c1=0, c2=0;
        for (double multi : multis)
        {
          Range<Double> r1 = Tools.bound(mBound,  pts[i] -mSpan*multi* step, pts[i] +mSpan*multi* step),
                        r2 = Tools.bound(mBoundX, ptsX[j]-mSpan*multi*xstep, ptsX[j]+mSpan*multi*xstep);
          c1=Tools.sliceCounts(mPositives, r1, r2);
          c2=Tools.sliceCounts(mNegatives, r1, r2);
          if (c1+c2>minN) break; // modest requirement to ensure stat validity
        }
        double p = Double.NaN;
        if (c1+c2>minN)
        {
          double pd_pos=c1/base_pos, pd_neg=c2/base_neg;

          pd_pos = pd_pos < 0 ? 0 : pd_pos;
          pd_neg = pd_neg < 0 ? 0 : pd_neg;

          p = (pd_pos == 0 && pd_neg == 0) ? 0d : pd_pos * p_pos / (pd_pos * p_pos + pd_neg * p_neg);
        }

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

    try
    {
      return mTransition.value(data,x);
    }
    catch (OutOfRangeException oe)
    {
      return Double.NaN;
    }
  }
  public StringBuffer df(double slice)
  {
    StringBuffer buf = new StringBuffer();

    if (mTransition==null) toTransition(0.5,0.5,false);

    buf.append("tag\tdata\tx\tval\n");

    double step = (mBound.upperEndpoint()-mBound.lowerEndpoint())/slice, xstep = (mBoundX.upperEndpoint()-mBoundX.lowerEndpoint())/slice;
    // dump the pos and neg
    buf = df(buf, mPositives, slice,step,xstep, "Y");
    buf = df(buf, mNegatives, slice,step,xstep, "N");

    // lookup the transition
    for (int i=0; i<slice; i++)
    {
      double data = mBound.lowerEndpoint() + i*step;
      for (int j=0; j<slice; j++)
      {
        double x = mBoundX.lowerEndpoint() + j*xstep;
        buf.append("t"+"\t"+data+"\t"+x+"\t"+lookup(data,x)+"\n");
      }
    }

    return buf;
  }
  private StringBuffer df(StringBuffer buf, TreeBasedTable<Double, Double, Integer> data, double slice, double step, double xstep, String tag)
  {
    for (int i=0; i<slice; i++)
    {
      double d0 = mBound.lowerEndpoint() + i*step;
      for (int j=0; j<slice; j++)
      {
        double x0 = mBoundX.lowerEndpoint() + j*xstep;
        buf.append(tag+"\t"+d0+"\t"+x0+"\t"+Tools.sliceCounts(data,Range.closed(d0,d0+step), Range.closed(x0,x0+xstep))+"\n");
      }
    }
    return buf;
  }
}
