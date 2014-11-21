package org.ms2ms.alg;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import com.google.common.collect.TreeMultimap;
import org.apache.hadoop.hbase.util.Hash;
import org.expasy.mzjava.core.ms.peaklist.PeakAnnotation;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.RetentionTime;
import org.expasy.mzjava.core.ms.spectrum.RetentionTimeList;
import org.ms2ms.utils.Tools;

import java.util.Collection;
import java.util.Iterator;

/**
 * Created with IntelliJ IDEA.
 * User: hliu
 * Date: 8/23/14
 * Time: 5:54 AM
 * To change this template use File | Settings | File Templates.
 */
public class Spectra
{
  public static boolean before(RetentionTimeList rt, double limit)
  {
    if (rt==null) return false;

    Iterator<RetentionTime> itr = rt.listIterator();
    while (itr.hasNext())
    {
      RetentionTime r = itr.next();
      // quit if any lower bound equals or exceeds the limit
      if (r.getMinRetentionTime()/60d>=limit) return false;
    }
    return true;
  }
  public static boolean after(RetentionTimeList rt, double limit)
  {
    if (rt==null) return false;

    Iterator<RetentionTime> itr = rt.listIterator();
    while (itr.hasNext())
    {
      RetentionTime r = itr.next();
      // quit if any upper bound equals or less than the limit
      if (r.getMaxRetentionTime()/60d<=limit) return false;
    }
    return true;
  }

  public static boolean contains(RetentionTimeList rt, Range<Double> bound)
  {
    if (rt==null || bound==null) return false;

    Iterator<RetentionTime> itr = rt.listIterator();
    while (itr.hasNext())
    {
      RetentionTime r = itr.next();
      if (bound.contains(r.getMaxRetentionTime()/60d) ||
          bound.contains(r.getMinRetentionTime()/60d)) return true;
    }
    return false;
  }
  public static <S extends PeakList<PeakAnnotation>> Multimap<Integer, S> toChargePeakList(Collection<S> spectra)
  {
    if (!Tools.isSet(spectra)) return null;

    Multimap<Integer, S> z_spec = HashMultimap.create();
    for (S spec : spectra) z_spec.put(spec.getPrecursor().getCharge(), spec);

    return z_spec;
  }
}
