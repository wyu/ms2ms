package org.ms2ms.data.ms;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/** A 'tall' array of sortable objects as a more space efficient collection vs SortedMultimap **/
/*  It's made of vertical piles of arrays. The horizontal is then called 'trunk' to mark the  */
/*  the beg and end of section in use. */
abstract public class AbstractPile<T extends Comparable<T>>
{
  protected int mInitialLength   = 100000, // the length of the array in each trunk
                mSeriesCapacity  = 255,
                mTrunkEnd        = 0,      // last of the pile in use
                mSeriesEnd       = 0;      // end position of the FragmentEntry of a key block

  protected Integer curr         = null,     // currnt trunk
                    currKey      = null;  // the key block we're looking at

  protected List<T[]> mDataPiles = new ArrayList<>();
  protected PileTrunk[]   trunks = new PileTrunk[mSeriesCapacity];
  protected T[]       mSeries;

  public AbstractPile()       { super(); }
  public AbstractPile(int s)  { super(); mInitialLength=s; }

  abstract public int    getKeyAt(int pile, int idx);
  abstract protected T[] newPile(int s);

  public int     getTrunkEnd()    { return mTrunkEnd; }
  public Integer getCurrentKey()  { return currKey; }
  public Integer getCurrPile()    { return curr; }

  public T get(int pile, int idx) { return mDataPiles.get(pile)[idx]; }
  public T at(int i)              { return mSeries[i]; }

  public T add(T s)
  {
    // do we enough room to populate it?
    hasNext();
    mDataPiles.get(curr)[trunks[curr].N++] = s;
//    mDataPiles.get(curr)[trunks[curr].N++] = new FragmentMatch(s, i);
    return s;
  }
  public T get(int idx) { return mDataPiles.get(curr)[idx]; }
  public T front()      { return mDataPiles.get(curr)[trunks[curr].beg]; }
  public T back()       { return mDataPiles.get(curr)[trunks[curr].N]; }

  public AbstractPile sort()
  {
    for (int i=0; i<getTrunkEnd(); i++)
      Arrays.parallelSort(mDataPiles.get(i), 0, trunks[i].N);

    return this;
  }
  public int size()
  {
    int n = 0;
    for (int i=0; i<getTrunkEnd(); i++) n += trunks[i].N;

    return n;
  }
  // prepare for the start of key iteration, after the sort
  public void start()
  {
    mSeriesEnd=0; currKey=0; curr=0;
    for (int i=0; i<getTrunkEnd(); i++)
      trunks[i].key = getKeyAt(i,trunks[i].beg);
  }
  // reset all the pointers back to an empty pile
  public AbstractPile init()
  {
    // need at least one trunk
    if (getTrunkEnd()==0) mDataPiles.add(newPile(mInitialLength));

    mSeries = newPile(mSeriesCapacity);
    mSeriesEnd=mTrunkEnd=curr=0;

    trunks[curr] = new PileTrunk();
    return this;
  }

  // make sure we have the place for the new comer
  public AbstractPile hasNext()
  {
    // do we enough room to populate it?
    if (trunks[curr].N>=mInitialLength) curr++;
    if (curr>=getTrunkEnd())
    {
      mDataPiles.add(newPile(mInitialLength));
      mTrunkEnd++;
      // reset the counts
      trunks[curr] = new PileTrunk();
    }

    return this;
  }

  public Integer nextStart()
  {
    Integer newStart=null;
    for (int i=0; i<getTrunkEnd(); i++)
      // the keys are sorted in desending order, so we want the largest one across the piles
      if (trunks[i].key!=null && (newStart==null || trunks[i].key>newStart))
        newStart=trunks[i].key;

    return newStart;
  }
  // return the chunk of FragmentEntry for the next key
  public int nextChunk()
  {
    currKey = nextStart(); mSeriesEnd=0;

//    System.out.print("Searching for "+currKey+" --> ");

    // quit if no more key to be find
    if (currKey==null) return 0;

    int i=0,j=0;
    for (i=0; i<getTrunkEnd(); i++)
    {
      // look for the entries for this peptides across the trunks
      boolean OK=false;
      if (getKeyAt(i,trunks[i].beg)==currKey)
      {
        for (j=trunks[i].beg; j<trunks[i].N; j++)
          if (getKeyAt(i,j)==currKey)
          {
            // accumualte the entries
            mSeries[mSeriesEnd++] = get(i,j);
          }
          else
          {
            trunks[i].beg =j;
            trunks[i].key =getKeyAt(i,j);
            OK=true; break;
          }
        // hit the end of the array without seeing a new peptide
        if (j>=trunks[i].N && !OK) trunks[i].key=null;
      }
    }
    return mSeriesEnd;
  }
}
