package org.ms2ms.data.collect;

import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;
import org.ms2ms.utils.Tools;

import java.util.*;

/**
 * ** Copyright 2014-2015 ms2ms.org
 * <p/>
 * Description:
 * <p/>
 * Author: wyu
 * Date:   11/1/14
 */
public class MapOfMultimap<K extends Comparable, L extends Comparable, T extends Comparable>
{
  private SortedMap<K, TreeMultimap<L, T>> mData;

  public MapOfMultimap() { super(); mData = new TreeMap<K, TreeMultimap<L, T>>(); }

  public long size()
  {
    if (!Tools.isSet(mData)) return 0;

    long counts = 0;
    for (Multimap<L,T> ts : mData.values()) counts += ts.size();
    return counts;
  }
//  public long key_size()
//  {
//    long counts = 0;
//    for (Multimap<L,T> ts : getData().values()) counts += ts.getData().size();
//
//    return counts;
//  }
  public MapOfMultimap<K,L,T>  put(   K key, L lable,      T  data) { initByKey(key).put(lable,    data); return this; }
  public MapOfMultimap<K,L,T>  putAll(K key,   Multimap<L, T> data) { initByKey(key).putAll(       data); return this; }
  public MapOfMultimap<K,L,T>  putAll(K key, L lable, List<T> data) { initByKey(key).putAll(lable, data); return this; }

//  /** Overwrite the existing entry
//   *
//   * @param key
//   * @param data
//   */
//  public void set(K key, L lable, T data)
//  {
//    if (mData == null) mData = new TreeMap<K, MultiTreeMap<L, T>>();
//    MultiTreeMap<L, T> D = mData.cells(key);
//    if (D == null)
//    {
//      D = new MultiTreeMap<L, T>();
//      mData.put(key, D);
//    }
//    else D.clear();
//    D.add(lable, data);
//  }
  public Collection<K> keySet() { return mData.keySet(); }
  public Collection<T> values()
  {
    Collection<T> made = new ArrayList<T>();
    for (Multimap<L,T> ts : mData.values())
      for (T t : ts.values()) made.add(t);
    return made;
  }
  public Collection<L> labelSet()
  {
    Collection<L> labels = new HashSet<L>();
    for (K k : keySet())
      labels.addAll(get(k).keySet());
    return labels;
  }

  public Multimap<L,T> get(K key)          { return mData==null ? null : mData.get(key); }
  public Collection<T> get(K key, L lable) { return mData==null || mData.get(key)==null ? null : mData.get(key).get(lable); }

  public SortedMap<K, TreeMultimap<L, T>> getData() { return mData; }
  public void clear() { if (getData() != null) getData().clear(); }

  public boolean remove(K key, L label, T val)
  {
    if (val == null) return false;

    boolean outcome = get(key, label).remove(val);
    if (outcome)
    {
      if (!Tools.isSet(get(key, label))) get(key).keySet().remove(label);
      if (!Tools.isSet(get(key)))        getData().remove(key);
    }
    return outcome;
  }
  public boolean remove(K key1, K key2)
  {
    Set<K> removed = new HashSet<K>(getData().subMap(key1, key2).keySet());

    if (Tools.isSet(removed))
    {
      for (K key : removed) getData().keySet().remove(key);
      return true;
    }
    removed = null;

    return false;
  }
  public Collection<T> subset(K k0, K k1, L l0, L l1)
  {
    Collection<TreeMultimap<L, T>> slice = mData.subMap(k0, k1).values();
    Collection<T>                 values = new HashSet<T>();
    if (Tools.isSet(slice))
      for (TreeMultimap<L, T> sub : slice)
      {
        Map<L, Collection<T>> s = sub.asMap().subMap(l0, l1);
        if (s!=null)
          for (Collection<T> t : s.values()) values.addAll(t);
      }

    return values;
  }
  public Collection<T> getByLabel(L s)
  {
    Collection<T> vals = new ArrayList<T>();
    for (K k : keySet())
      if (get(k, s) != null)
        vals.addAll(get(k, s));

    return vals;
  }

  /** make sure it's ready to accept the insertion by key
   *
   * @param key
   * @return
   */
  private Multimap<L, T> initByKey(K key)
  {
    if (mData == null) mData = new TreeMap<K, TreeMultimap<L, T>>();
    TreeMultimap<L, T> D = mData.get(key);
    if (D == null) {
      D = TreeMultimap.create();
      mData.put(key, D);
    }
    return D;
  }
  public static <R extends Comparable, C extends Comparable, V extends Comparable> MapOfMultimap<R,C,V>
  create() { return new MapOfMultimap<>(); }
}
