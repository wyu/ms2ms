package org.ms2ms.utils;

import com.google.common.collect.Range;

import java.util.*;

/**
 * Created by wyu on 4/17/14.
 */
public class Tools
{
  static String ZEROES = "000000000000";
  static String BLANKS = "            ";

  public static     boolean isSet(String        s) { return s!=null && s.length()>0; }
  public static <T> boolean isSet(Collection<T> s) { return s!=null && s.size()>0; }
  public static <T> boolean isSet(T[]           s) { return s!=null && s.length>0; }
  public static     boolean isSet(int[]         s) { return s!=null && s.length>0; }

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

  //	Author:  Jean Vaucher
  //  Date:    March 2006
  // http://www.iro.umontreal.ca/~vaucher/Java/tutorials/Formatting.html
  /* -------------------------------------------------------------------------
     Meant to be used in a "print" statement to align data (text, int or real)
     in fields of W characters: on the right if W>0 and left if W<0.
     With real numbers, the parameter N is the number of decimals required.

         Util.format( String data, W )
         Util.format( int data, W )
         Util.format( double data, N, W )

     Usage:
         System.out.println("Name            age  result");
         System.out.println( Util.format(name,-15)
           + Util.format(age,5)
           + Util.format(mark,2,7));
     ------------------------------------------------------------------------

     - In my "formatting" page, "format" used to be called "toString"

     ------------------------------------------------------------------------ */
  public static String format( int val, int w) {	return format( Integer.toString(val), w); }
  public static String format( String s, int w)
  {
    int w1 = Math.abs(w);
    int n = w1-s.length();

    if ( n <= 0 ) return s;
    while (BLANKS.length()<n) BLANKS += "      ";
    if ( w<0 )
      return s + BLANKS.substring(0,n);
    else
      return BLANKS.substring(0,n) + s;
  }

  public static String format( double val, int n, int w)
  {
    //	rounding
    double incr = 0.5;
    for( int j=n; j>0; j--) incr /= 10;
    val += incr;

    String s = Double.toString(val);
    int n1 = s.indexOf('.');
    int n2 = s.length() - n1 - 1;

    if (n>n2)  {
      int len = n-n2;
      while (ZEROES.length()<len) ZEROES += "000000";
      s = s+ZEROES.substring(0, len);
    }
    else if (n2>n) s = s.substring(0,n1+n+1);

    return format( s, w );
  }
  public static String d2s(double s, int i)
  {
    return format(s, i, 0);
  }

  public static <C extends Comparable> Collection<Range<C>> merge(Collection<Range<C>> r)
  {
    // merge the slices
    Collection<Range<C>> pool = new ArrayList<Range<C>>();
    Set<Range<C>>   discarded = new HashSet<Range<C>>();
    for (Range<C> r1 : r)
      if (!discarded.contains(r1))
      {
        Range<C> p = r1;
        for (Range<C> r2 : r)
          if (r2 != r1 && !discarded.contains(r2) && r1.isConnected(r2)) { p = p.span(r2); discarded.add(r2); }
        pool.add(p);
      }

    return pool;
  }
  public static String extend(String s0, String s1, String delimiter)
  {
    return s0==null?s1:(s0+delimiter+s1);
  }
  public static <K, V> Map<K,V> putAll(Map<K,V> map, Map<K,V> in)
  {
    if (map!=null && in!=null) map.putAll(in);
    return map;
  }
  public static <T> boolean has(T[] vals, T obj, int end)
  {
    if (isSet(vals) && obj!=null)
    {
      int bound = (end>=0&&end<=vals.length?end:vals.length);
      for (int i=0; i<bound; i++)
        if (vals[i]!=null && vals[i].equals(obj)) return true;
    }

    return false;
  }
}
