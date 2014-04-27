package org.ms2ms.utils;

import java.util.Collection;

/**
 * Created by wyu on 4/26/14.
 */
public class Stats
{
  public static double mean(Collection<Double> s)
  {
    if (!Tools.isSet(s)) return 0;

    double avg = 0d;
    for (Double v : s) avg+=v;
    return avg/(double )s.size();
  }
  public static double stdev(Collection<Double> s) { return Math.sqrt(variance(s)); }
  public static double variance(Collection<Double> s)
  {
    if (!Tools.isSet(s) || s.size()==1) return 0;

    double avg=mean(s), var=0d;
    for (Double v : s) { var+= (v-avg)*(v-avg); }

    return var/((double )s.size()-1d);
  }
}
