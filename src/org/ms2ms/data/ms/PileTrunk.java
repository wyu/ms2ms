package org.ms2ms.data.ms;

// keeping the stats on each of the trunk
public class PileTrunk
{
  public int N,             // total number of actual element in this trunk
             beg;           // beginning and end of the key block we're looking at
  public int key=Integer.MAX_VALUE;  // the key index of the block we're looking at

  public PileTrunk() { N=beg=0; }

  public PileTrunk setKey(int s) { key=s; return this; }
}
