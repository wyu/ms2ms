package org.ms2ms.alg;

import com.google.common.collect.Range;
import org.ms2ms.data.NameValue;
import org.ms2ms.utils.Stats;
import org.ms2ms.utils.Tools;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Comments:
 */
public class Parser
{
/** Parse comma-separated values (CSV), a common Windows file format.
 * Sample input: "LU",86.25,"11/4/1998","2:19PM",+4.0625
 * <p>
 * Inner logic adapted from a C++ original that was
 * Copyright (C) 1999 Lucent Technologies
 * Excerpted from 'The Practice of Programming'
 * by Brian W. Kernighan and Rob Pike.
 * <p>
 * Included by permission of the http://tpop.awl.com/ web site,
 * which says:
 * "You may use this code for any purpose, as long as you leave
 * the copyright notice and book citation attached." I have done so.
 * author: Brian W. Kernighan and Rob Pike (C++ original)
 * author: Ian F. Darwin (translation into Java and removal of I/O)
 * author: Ben Ballard (rewrote advQuoted to handle '""' and for readability)
 */

  public static final char COMMA = ',';
  public static final char TABB  = '\t';
  public static final char SPACE = ' ';
  public static final char QUOTE2 = '"';

  /** parse: break the input String into fields
   * @return java.util.Iterator containing each field
   * from the original as a String, in order.
   */
  public static List<String> parse(String line, char fieldSep)
  {
    StringBuffer sb   = new StringBuffer();
    List<String> list = new ArrayList<String>();

    int i = 0;
    if (line.length() == 0)
    {
      list.add(line);
      return list;
    }

    do
    {
      sb.setLength(0);
      if (i < line.length() && line.charAt(i) == QUOTE2) i = advQuoted(line, sb, ++i, fieldSep);	// skip quote
      else                                               i = advPlain(line, sb,    i, fieldSep);
      list.add(sb.toString());
      i++;
    } while (i < line.length());

    return list;
  }

  /** advQuoted: quoted field; return index of next separator */
  protected static int advQuoted(String s, StringBuffer sb, int i, char fieldSep)
  {
    int j;
    int len= s.length();
    for (j=i; j<len; j++)
    {
      if (s.charAt(j) == QUOTE2 && j+1 < len)
      {
        if (s.charAt(j+1) == QUOTE2)
        {
          j++; // skip escape char
        }
        else if (s.charAt(j+1) == fieldSep)
        { //next delimeter
          j++; // skip end quotes
          break;
        }
      }
      else if (s.charAt(j) == QUOTE2 && j+1 == len)
      { // end quotes at end of line
        break; //done
      }
      sb.append(s.charAt(j));	// regular character.
    }
    return j;
  }

  /** advPlain: unquoted field; return index of next separator */
  protected static int advPlain(String s, StringBuffer sb, int i, char fieldSep)
  {
    int j = s.indexOf(fieldSep, i); // look for separator
    int q = s.indexOf(QUOTE2, i); // look for separator

    if (q != -1 && q < j)
    {
      q = s.indexOf(QUOTE2, q+1); // move to the next one
      j = s.indexOf(fieldSep, q); // look for separator
    }

    //log.debug("csv: " + "i = " + i + " j = " + j);
    if (j == -1)
    {               	// none found
      sb.append(s.substring(i));
      return s.length();
    } else {
      sb.append(s.substring(i, j));
      return j;
    }
  }
  public static Map<String, String> extractVariables(String pattern, File[] files)
  {
    String left = null, right = null;

    if (Tools.isSet(pattern))
    {
      left  = pattern.substring(0, pattern.indexOf("*"));
      right = pattern.substring(pattern.lastIndexOf("*") + 1, pattern.length());
    }

    Map<String, String> outs = new HashMap<String, String>();
    for (File F : files)
    {
      //if (!F.getPath().matches(pattern)) continue;
      String L = F.getPath();
      if (Tools.isSet(pattern))
      {
        int i_left = L.indexOf(left), i_right = L.lastIndexOf(right);
        if (i_left == 0 && i_right + right.length() == L.length())
          outs.put(L.substring(i_left + left.length(), i_right), L);
      }
      else outs.put(L, L);
    }
    return outs;
  }
  public static NameValue newNameValue(String line, char fieldSep)
  {
    if (line == null) return null;
    
    NameValue nv = new NameValue();
    int      pos = line.indexOf(fieldSep);

    if (pos < 0) { nv.name = line; return nv; }

    nv.name = line.substring(0,  pos).trim();
    nv.val  = line.substring(pos + 1).trim();
    return nv;
  }
  public static NameValue newNameValue(String line, String fieldSep)
  {
    if (line == null) return null;

    NameValue nv = new NameValue();
    int      pos = line.indexOf(fieldSep);

    if (pos < 0) { nv.name = line; return nv; }

    nv.name = line.substring(0,  pos                ).trim();
    nv.val  = line.substring(pos + fieldSep.length()).trim();
    return nv;
  }
  public static Range<Double> newRange(String line, String fieldSep)
  {
    if (line == null) return null;

    int pos = line.indexOf(fieldSep);
    if (pos < 0) { return null; }

    return Range.closed(Stats.toDouble(line.substring(0,  pos                ).trim()),
                        Stats.toDouble(line.substring(pos + fieldSep.length()).trim()));
  }
}
