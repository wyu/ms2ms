package org.ms2ms.data.ms;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DoublePile
{
  protected int mInitialLength   = 100000, // the length of the array in each trunk
      mSeriesCapacity  = 255,
      mTrunkEnd        = 0;      // last of the pile in use

  protected Integer curr         = null,     // currnt trunk
      currKey      = null;  // the key block we're looking at

  protected List<double[]> mDataPiles = new ArrayList<>();
  protected PileTrunk[]   trunks = new PileTrunk[mSeriesCapacity];

  public DoublePile()       { super(); }
  public DoublePile(int s)  { super(); mInitialLength=s; init(); }

  public int     getTrunkEnd()    { return mTrunkEnd; }
  public Integer getCurrentKey()  { return currKey; }
  public Integer getCurrPile()    { return curr; }

  public double get(int pile, int idx) { return mDataPiles.get(pile)[idx]; }

  public DoublePile add(double s)
  {
    // do we enough room to populate it?
    hasNext();
    mDataPiles.get(curr)[trunks[curr].N++] = s;
//    mDataPiles.get(curr)[trunks[curr].N++] = new FragmentMatch(s, i);
    return this;
  }
  public double get(int idx) { return mDataPiles.get(curr)[idx]; }
  public double front()      { return mDataPiles.get(curr)[trunks[curr].beg]; }
  public double back()       { return mDataPiles.get(curr)[trunks[curr].N]; }

  public DoublePile sort()
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
    currKey=0; curr=0;
    for (int i=0; i<getTrunkEnd(); i++)
      trunks[i].key = 0;
  }
  // reset all the pointers back to an empty pile
  public DoublePile init()
  {
    // need at least one trunk
    if (getTrunkEnd()==0) mDataPiles.add(new double[mInitialLength]);

    mTrunkEnd=curr=0;

    trunks[curr] = new PileTrunk();
    return this;
  }

  // make sure we have the place for the new comer
  public DoublePile hasNext()
  {
    // do we enough room to populate it?
    if (trunks[curr].N>=mInitialLength) curr++;
    if (curr>=getTrunkEnd())
    {
      mDataPiles.add(new double[mInitialLength]);
      mTrunkEnd++;
      // reset the counts
      trunks[curr] = new PileTrunk();
    }

    return this;
  }
}
