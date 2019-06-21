package org.ms2ms.data.ms;

import com.google.common.collect.Ordering;
import org.ms2ms.data.Point;
import org.ms2ms.math.Points;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Tools;

import java.util.ArrayList;
import java.util.List;

public class HistogramPile extends DoublePile
{
  private Double mCentroid, mPct1, mPct5, mPct33, mMedian;
  private List<Point> mHistogram;

  public HistogramPile()       { super(); }
  public HistogramPile(int s)  { super(s); }

  public Double getCentroid()  { return mCentroid; }
  public Double getMedian()    { return mMedian; }
  public Double getTop1Pct()   { return mPct1; }
  public Double getTop5Pct()   { return mPct5; }
  public Double getTop33Pct()  { return mPct33; }

  // calc the distribution and stats
  public HistogramPile survey()
  {
    if (size()==0) return this;

    sort();

    double[][] stats = new double[2][getTrunkEnd()];
    List<Double> tops = new ArrayList<>();

    for (int t=0; t<getTrunkEnd(); t++)
    {
      // gather the stat in each pile first
      stats[0][t] = get(t, (int )Math.floor(0.50d*trunks[t].N));
      stats[1][t] = get(t, (int )Math.floor(0.67d*trunks[t].N));

      int cut = (int )Math.floor(0.95d*trunks[t].N);
      for (int i=cut; i<trunks[t].N; i++) tops.add(get(t,i));
    }

    // grab the percentiles
    mMedian = Stats.mean0(stats[0]);
    mPct33  = Stats.mean0(stats[1]);

    // sort the top ranked
    tops.sort(Ordering.natural().reversed());

    // look into the top ranks values for true top-percentiles
    mPct5   = tops.get((int )Math.floor((size()-1)*0.05d));
    mPct1   = tops.get((int )Math.floor((size()-1)*0.01d));

    tops = (List )Tools.dispose(tops);

    if (Tools.isSet(mHistogram))
    {
      mCentroid = Points.centroid(mHistogram);
    }
    return this;
  }



  // generate the histogram
  public HistogramPile generate(int step_num)
  {
    if (size()==0 || step_num == 0) return this;

    mHistogram = new ArrayList<>();

    sort();

    // either equal space or numbers of data points
    double N = size(), counts=0, sum=0d;

    for (int t=0; t<getTrunkEnd(); t++)
    {
//      int span = Math.round(trunks[t].N/step_num);
//      for (int j=0; j<trunks[t].N; j++)
//      {
//        sum+=get(t,j);
//        if (++counts>=span)
//        {
//          mHistogram.add(new Point(sum/counts, counts));
//          counts=0; sum=0d;
//        }
//      }
    }

    return this;
  }

}
