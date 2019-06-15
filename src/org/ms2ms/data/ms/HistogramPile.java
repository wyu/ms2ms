package org.ms2ms.data.ms;

import com.google.common.collect.Range;
import org.ms2ms.data.Point;
import org.ms2ms.math.Histogram;
import org.ms2ms.utils.Tools;

import java.util.Collections;
import java.util.List;

public class HistogramPile extends DoublePile
{
  private List<Point> mHistogram;

  public HistogramPile()       { super(); }
  public HistogramPile(int s)  { super(s); }

  // calc the distribution and stats




  // generate the histogram
  public HistogramPile generate(int step_num)
  {
    if (size()==0 || step_num == 0) return this;

    return this;
  }

}
