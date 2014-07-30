package org.ms2ms.utils;

import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 7/19/14
 * Time: 6:20 PM
 * To change this template use File | Settings | File Templates.
 */
public class ArrayKey
{
  public Object[] key;

  public ArrayKey(Object[] s) { super(); key=s; }

  @Override
  public boolean equals(Object obj)
  {
    if(obj!=null && obj instanceof ArrayKey)
    {
      ArrayKey s = (ArrayKey)obj;
      return Arrays.equals(key, s.key);
    }
    return false;
  }

  @Override
  public int hashCode()
  {
    return Strs.toString(key, "").hashCode();
  }
}
