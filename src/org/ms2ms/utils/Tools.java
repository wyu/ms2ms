package org.ms2ms.utils;

import java.util.Collection;
import java.util.Iterator;

/**
 * Created by wyu on 4/17/14.
 */
public class Tools
{
  public static <T> boolean isSet(Collection<T> s) { return s!=null && s.size()>0; }
  public static <T> T front(Collection<T> s)
  {
    if (isSet(s))
      for (T t : s) return t;

    return null;
  }
  public static <T> T back(Collection<T> s)
  {
    T last = null;
    if (isSet(s))
      for (T t : s) last=t;

    return last;
  }
}
