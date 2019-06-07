package org.ms2ms.data.ms;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/** A 'tall' array of sortable objects as a more space efficient collection vs SortedMultimap **/
//
abstract public class AbstractPile<T extends Comparable<T>>
{
  // keeping the stats on each of the trunk
  class Trunk
  {
    public int N,          // total number of actual element in this trunk
        beg, end,   // beginning and end of the key block we're looking at
        key;        // the key index of the block we're looking at

    public Trunk() { N=beg=end= key =0; }
  }

  protected int mInitialLength  = 100000, // the length of the array in each trunk
                mSeriesCapacity = 255,
                curr=0,                  // currnt trunk
                currKey=0,               // the key block we're looking at
                mSeriesEnd=0;            // end position of the FragmentEntry of a key block

  protected List<T[]> mDataPiles = new ArrayList<>();
  protected Trunk[]   trunks     = (Trunk[] )new AbstractPile.Trunk[mSeriesCapacity];
  protected T[]       mSeries;

  public AbstractPile()       { super(); }
  public AbstractPile(int s)  { super(); mInitialLength=s; }

  abstract public int    getKeyAt(int pile, int idx);
  abstract protected T[] newPile();

  public T get(int pile, int idx) { return mDataPiles.get(pile)[idx]; }

  public AbstractPile add(T s)
  {
    // do we enough room to populate it?
    hasNext();
    mDataPiles.get(curr)[trunks[curr].N++] = s;
//    mDataPiles.get(curr)[trunks[curr].N++] = new FragmentMatch(s, i);
    return this;
  }
  public T get(int idx) { return mDataPiles.get(curr)[idx]; }
  public T front()      { return mDataPiles.get(curr)[trunks[curr].beg]; }
  public T back()       { return mDataPiles.get(curr)[trunks[curr].end]; }

  public AbstractPile sort()
  {
    for (int i=0; i<mDataPiles.size(); i++)
      Arrays.parallelSort(mDataPiles.get(i), 0, trunks[i].N);

    return this;
  }
  public int size()
  {
    int n = 0;
    for (int i=0; i<mDataPiles.size(); i++) n += trunks[i].N;

    return n;
  }
  // prepare for the start of key iteration, after the sort
  public void start()
  {
    mSeriesEnd=0; currKey=0;
    for (int i=0; i<mDataPiles.size(); i++)
      trunks[i].key = getKeyAt(i,trunks[i].beg);
  }
  // reset all the pointers back to an empty pile
  public AbstractPile init()
  {
    // need at least one trunk
    if (mDataPiles.size()==0) mDataPiles.add(newPile());

    mSeriesEnd=curr=0;
    Arrays.fill(trunks, new Trunk());
    return this;
  }

  // make sure we have the place for the new comer
  public AbstractPile hasNext()
  {
    // do we enough room to populate it?
    if (trunks[curr].end>mInitialLength) curr++;
    if (curr>=mDataPiles.size())
    {
      mDataPiles.add(newPile());
      curr++; // move to the new trunk
    }

    return this;
  }

  public int nextStart(int currStart)
  {
    int newStart=currStart;
    for (int i=0; i<mDataPiles.size(); i++)
      if (newStart==0 || trunks[i].key >currKey) newStart=trunks[i].key;

    return newStart;
  }
  // return the chunk of FragmentEntry for the next key
  public int nextPeptide()
  {
    mSeriesEnd=0; currKey = nextStart(currKey);
    for (int i=0; i<mDataPiles.size(); i++)
    {
      // look for the entries for this peptides across the trunks
      if (getKeyAt(i,trunks[i].beg)==currKey)
        for (int j=trunks[i].beg; j<trunks[i].N; j++)
          if (getKeyAt(i,j)==currKey)
          {
            // accumualte the entries
            mSeries[mSeriesEnd++] = get(i,j);
          }
          else
          {
            trunks[i].beg=j; trunks[i].key =getKeyAt(i,j);
            break;
          }
    }
    return mSeriesEnd;
  }


}
