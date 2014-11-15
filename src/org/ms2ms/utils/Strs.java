package org.ms2ms.utils;

import com.google.common.collect.Table;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 6/13/14
 * Time: 8:27 AM
 * To change this template use File | Settings | File Templates.
 */
public class Strs
{
  public static String extend(String s0, String s1, String delimiter)
  {
    return s0==null?s1:(s0+delimiter+s1);
  }
  public static String dequotes(String s, char q)
  {
    if (s.length() > 0              &&
        s.charAt(0)            == q &&
        s.charAt(s.length()-1) == q)
      return s.length() - 1 > 1 ? s.substring(1, s.length() - 1) : s.substring(1);
    return s;
  }
  public static void dequotes(String[] ss, char q)
  {
    for (int i = 0; i != ss.length; i++)
      ss[i] = dequotes(ss[i], q);
  }
  public static String[] trim(String[] ss)
  {
    for (int i = 0; i != ss.length; i++)
      ss[i] = new String(ss[i].trim());

    return ss;
  }

  public static String trim(String s) { return s!=null?s.trim():s; }
  public static List<String> split(String s, String regx, boolean trim)
  {
    List<String> list  = new ArrayList<String>();
    String[]     items = s.split(regx);
    for (String ss : items)
      if (!trim || Tools.isSet(ss.trim()))
        list.add(trim ? ss.trim() : ss);
    // return the list
    return list;
  }
  public static String[] split(String s, char c, boolean trim)
  {
    int i, b, e;
    int cnt;
    String res[];
    int ln = s.length();

    i = 0;
    cnt = 1;
    while ((i = s.indexOf(c, i)) != -1) {
      cnt++;
      i++;
    }
    res = new String[cnt];

    i = 0;
    b = 0;
    while (b <= ln) {
      e = s.indexOf(c, b);
      if (e == -1) {
        e = ln;
      }
      if (!trim) {
        res[i++] = s.substring(b, e);
      } else {
        int e2 = e - 1;
        while (e2 >= 0 && Character.isWhitespace(s.charAt(e2))) {
          e2--;
        }
        e2++;
        while (b < ln && Character.isWhitespace(s.charAt(b))) {
          b++;
        }
        if (b < e) {
          res[i++] = s.substring(b, e2);
        } else {
          res[i++] = "";
        }
      }
      b = e + 1;
    }
    return res;
  }
  public static String[] split(String s, char t)
  {
    if (s == null || t == 0) return null;

    return split(s, t, false);
  }
  public static String toString(String[] ss, String dl) {
    return toString(ss, dl, 0, ss != null ? ss.length : 0);
    /*// no point to go further
    if (ss == null) return null;
    String made = new String();
    for (int i = 0; i != ss.length; i++) if (ss[i] != null) made = extend(made, ss[i], dl);
    // return the concatenated string
    return made; */
  }
  public static <K, T> String toString(Map<K, T> ss, String dl) {
    // no point to go further
    if (ss == null) return null;
    String made = new String();
    for (Map.Entry E : ss.entrySet())
      made = extend(made, E.getKey().toString() + "--" + E.getValue().toString(), dl);
    // return the concatenated string
    return made;
  }
  public static <T> String toString(T[] ss, String dl, Integer start, Integer end) {
    // no point to go further
    if (ss == null) return null;
    String made = new String();
    for (int i = start; i != Math.min(end, ss.length); i++)
      if (ss[i] != null) made = extend(made, ss[i].toString(), dl);
    // return the concatenated string
    return made;
  }
  public static <T> String toString(T[] ss, String dl) {
    // no point to go further
    if (ss == null) return null;
    String made = new String();
    for (int i = 0; i != ss.length; i++)
      if (ss[i] != null) made = extend(made, ss[i].toString(), dl);
    // return the concatenated string
    return made;
  }
  public static String toString(double[] ss, String dl) {
    // no point to go further
    if (ss == null) return null;
    String made = new String();
    for (int i = 0; i != ss.length; i++)
      made = extend(made, ss[i] + "", dl);
    // return the concatenated string
    return made;
  }
  public static <T> String toString(Collection<T> ss, String dl) {
    // no point to go further
    if (ss == null) return null;
    String made = new String();
    for (T t : ss)
      if (t != null) made = extend(made, t.toString(), dl);
    // return the concatenated string
    return made;
  }
  public static String toString(Collection<Double> ss, int deci, String dl) {
    // no point to go further
    if (ss == null) return null;
    String made = new String();
    for (Double t : ss)
      if (t != null) made = extend(made, Tools.d2s(t, deci), dl);
    // return the concatenated string
    return made;
  }
  public static String toString(Double[] ss, int deci, String dl) {
    // no point to go further
    if (ss == null) return null;
    String made = new String();
    for (Double t : ss)
      if (t != null) made = extend(made, Tools.d2s(t, deci), dl);
    // return the concatenated string
    return made;
  }
  public static String toString(String[] ss, String dl, Integer start, Integer end) {
    // no point to go further
    if (ss == null) return null;
    String made = new String();
    for (int i = start; i != end; i++) if (ss[i] != null) made = extend(made, ss[i], dl);
    // return the concatenated string
    return made;
  }
  public static <T> List<String> toStrings(List<T> ss) {
    // no point to go further
    if (ss == null) return null;
    List<String> made = new ArrayList<String>();
    for (T t : ss)
      if (t != null) made.add(t.toString());
    // return the concatenated string
    return made;
  }
/*
  public static String[] toStringArray(Collection<String> ss)
  {
    if (!Tools.isSet(ss)) return null;

    String[] cols = new String[ss.size()];
    int order = 0;
    for (String s : ss) cols[order++] = s;

    return cols;
  }
*/
  public static <T> Collection<String> toStrings(Collection<T> ss) {
    // no point to go further
    if (ss == null) return null;
    Collection<String> made = new ArrayList<String>();
    for (T t : ss)
      if (t != null) made.add(t.toString());
    // return the concatenated string
    return made;
  }
  public static String toQuotedString(String[] ss, String dl, String quote) {
    // no point to go further
    if (ss == null) return null;
    String made = new String();
    for (int i = 0; i != ss.length; i++)
      if (ss[i] != null) made = extend(made, quote + ss[i] + quote, dl);
    // return the concatenated string
    return made;
  }
  public static String toQuotedString(Collection<String> ss, String dl, String quote) {
    // no point to go further
    if (ss == null) return null;
    String made = new String();
    for (String s : ss)
      if (s != null) made = extend(made, quote + s + quote, dl);
    // return the concatenated string
    return made;
  }
  public static String toQuotedString(Object[] ss, String dl, String quote) {
    // no point to go further
    if (ss == null) return null;
    String made = new String();
    for (int i = 0; i != ss.length; i++)
      if (ss[i] != null) made = extend(made, quote + ss[i].toString() + quote, dl);
    // return the concatenated string
    return made;
  }
  public static String toString(List<String> ss, String dl) {
    // no point to go further
    if (ss == null) return null;
    String made = new String();
    for (String s : ss) made = extend(made, s, dl);
    // return the concatenated string
    return made;
  }
  public static String toStringFromLong(Collection<Long> ss, String dl) {
    // no point to go further
    if (ss == null) return null;
    String made = new String();
    for (Long s : ss) made = extend(made, s.toString(), dl);
    // return the concatenated string
    return made;
  }
  public static String toStringFromLong(Long[] ss, String dl) {
    // no point to go further
    if (ss == null) return null;
    String made = new String();
    for (Long s : ss) made = extend(made, s.toString(), dl);
    // return the concatenated string
    return made;
  }
  public static String[] toStringArray(Collection s)
  {
    String[] out = new String[s.size()];
    int i=0;
    for (Object v : s)  out[i++]=v.toString();

    return out;
  }
  public static String[] toStringArray(Object[] s)
  {
    String[] out = new String[s.length];
    int i=0;
    for (Object v : s)  out[i++]=v.toString();

    return out;
  }
  public static String toString(Table t)
  {
    if (!Tools.isSet(t)) return null;

    StringBuffer buf = new StringBuffer();
    for (Object col : t.columnKeySet()) buf.append(col.toString() + "\t");
    buf.append("\n");
    for (Object row : t.rowKeySet())
    {
      for (Object col : t.columnKeySet())
      {
        Object val = t.get(row, col);
        buf.append((val!=null?val.toString():"--") + "\t");
      }
      buf.append("\n");
    }
    return buf.toString();
  }
  public static String rtuncate(String s, int n)
  {
    return s==null||s.length()<n?s:s.substring(s.length()-n);
  }
  /** Converts time in milliseconds to a <code>String</code> in the format HH:mm:ss.SSS.
   *  lifted from http://www.uk-dave.com/bytes/java/long2time.php
   *
   * @param time the time in milliseconds.
   * @return a <code>String</code> representing the time in the format HH:mm:ss.SSS.
  */
  public static String msecToString(long time)
  {
    int milliseconds = (int)( time          % 1000);
    int      seconds = (int)((time/1000   ) % 60);
    int      minutes = (int)((time/60000  ) % 60);
    int        hours = (int)((time/3600000) % 24);
    String millisecondsStr = (milliseconds<10 ? "00" : (milliseconds<100 ? "0" : ""))+milliseconds;
    String      secondsStr = (seconds<10 ? "0" : "")+seconds;
    String      minutesStr = (minutes<10 ? "0" : "")+minutes;
    String        hoursStr = (hours<10 ? "0" : "")+hours;
    return new String(hoursStr+":"+minutesStr+":"+secondsStr+"."+millisecondsStr);
  }

}
