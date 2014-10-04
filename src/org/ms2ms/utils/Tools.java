package org.ms2ms.utils;

import com.google.common.collect.*;
import org.apache.commons.lang.NumberUtils;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.ms2ms.data.Dataframe;
import org.ms2ms.data.MultiTreeTable;
import org.ms2ms.data.Var;

import java.util.*;

/**
 * Created by wyu on 4/17/14.
 */
public class Tools
{
  static String ZEROES = "000000000000";
  static String BLANKS = "            ";

  public static     boolean isSet(String         s) { return s!=null && s.length()>0; }
  public static <T> boolean isSet(Collection<T>  s) { return s!=null && s.size()>0; }
  public static <T> boolean isSet(Map            s) { return s!=null && s.size()>0; }
  public static <T> boolean isSet(MultiTreeTable s) { return s!=null && s.size()>0; }
  public static <T> boolean isSet(TreeMultimap   s) { return s!=null && s.size()>0; }
  public static <T> boolean isSet(T[]            s) { return s!=null && s.length>0; }
  public static <T> boolean isSet(double[]       s) { return s!=null && s.length>0; }
  public static     boolean isSet(int[]          s) { return s!=null && s.length>0; }
  public static     boolean isSet(Table          s) { return s!=null && s.size()>0; }

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
  public static <T> T  front(     T... s) { return s!=null?s[0]:null; }
  public static <T> T  back(      T... s) { return s!=null?s[s.length-1]:null; }
  public static double front(double... s) { return s!=null?s[0]:null; }
  public static double back( double... s) { return s!=null?s[s.length-1]:null; }

  //	Author:  Jean Vaucher
  //  Date:    March 2006
  // http://www.iro.umontreal.ca/~vaucher/Java/tutorials/Formatting.html
  /* -------------------------------------------------------------------------
     Meant to be used in a "print" statement to run data (text, int or real)
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
  public static <K, V> Map<K,V> putAll(Map<K,V> map, Map<K,V> in)
  {
    if (map!=null && in!=null) map.putAll(in);
    return map;
  }
  public static boolean contains(int[] tt, int t)
  {
    if (isSet(tt))
      for (int i : tt)
        if (i==t) return true;

    return false;
  }
  public static <T> boolean contains(T[] vals, T obj, int end)
  {
    if (isSet(vals) && obj!=null)
    {
      int bound = (end>=0&&end<=vals.length?end:vals.length);
      for (int i=0; i<bound; i++)
        if (vals[i]!=null && vals[i].equals(obj)) return true;
    }

    return false;
  }
  public static <T extends Object> boolean isA(T A, T... B)
  {
    if ((A == null && B != null) ||
        (A != null && B == null)) return false;
    if (A == null && B == null) return true;

    for (T b : B) if (A.equals(b)) return true;
    return false;
  }
  public static Collection dispose(Collection s)
  {
    if (s!=null) { s.clear(); s=null; }
    return null;
  }
  public static Collection dispose(Map s)
  {
    if (s!=null) { s.clear(); s=null; }
    return null;
  }
  public static Collection dispose(Multimap s)
  {
    if (s!=null) { s.clear(); s=null; }
    return null;
  }
  public static Collection dispose(Table s)
  {
    if (s!=null) { s.clear(); s=null; }
    return null;
  }
  public static Map putNotNull(Map m, Object k, Object v)
  {
    if (m!=null && k!=null && v!=null) m.put(k, v);
    return m;
  }
  public static double[][] sort(double[]... xs)
  {
    // bubble sort from
    // http://thilinasameera.wordpress.com/2011/06/01/sorting-algorithms-sample-codes-on-java-c-and-matlab/
    int lenD = xs[0].length; double tmp = 0;
    for(int i = 0;i<lenD;i++)
      for(int j = (lenD-1);j>=(i+1);j--)
        if(xs[0][j]<xs[0][j-1])
          for (int k=0; k<xs.length; k++)
          {
            tmp = xs[k][j];
            xs[k][j]=xs[k][j-1];
            xs[k][j-1]=tmp;
          }

    return xs;
  }
  public static <T extends Object> Collection<T> slice(TreeBasedTable<Double, Double, T> data, Double r1, Double r2, Double c1, Double c2)
  {
    Collection<T> results = new ArrayList<T>();

    if (isSet(data) && (r2>r1) && (c2>c1))
    {
      SortedMap<Double,Map<Double, T>> s1 = data.rowMap().subMap(r1, r2);
      if (isSet(s1))
      {
        SortedMap<Double,Map<Double, T>> s2 = s1.subMap(c1, c2);
        if (isSet(s2))
          for (Map<Double, T> s3 : s2.values())
            results.addAll(s3.values());
        s2=null;
      }
      s1=null;
    }
    return results;
  }
  public static <K, V> BiMap<K, V> putNew(BiMap<K, V> m, K key, V val)
  {
    if (m==null || key==null || val==null) return m;
    // Check if the row or the aligned row is already filled
    if (!m.containsKey(key) && !m.inverse().containsKey(val)) m.put(key, val);

    return m;
  }
  public static <T> Collection<T> common(Collection<T> A, Collection<T> B)
  {
    if (A==null || B==null) return A;

    Iterator<T> itr = A.iterator();
    while (itr.hasNext())
      if (!B.contains(itr.next())) itr.remove();

    return A;
  }
  public static boolean equals(Object A, Object B) { return A!=null&&B!=null?A.equals(B):false; }
  public static boolean equalsCaseless(String A, String B) { return A!=null&&B!=null?A.equalsIgnoreCase(B):false; }
}
