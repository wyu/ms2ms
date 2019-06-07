package org.ms2ms.data.ms;

import java.lang.reflect.Array;

// TODO not completed yet!!
public class FixedArray<T>
{
  private T[] data;
  private int capacity=255, end=0;

  public FixedArray()      { }
  public FixedArray(int s, Class<T> C) { capacity=s; init(C); }

  private void init(Class<T> C) { data = (T[] )Array.newInstance(C, capacity); }

  public FixedArray add(T s) { data[end++]=s; return this; }

}
