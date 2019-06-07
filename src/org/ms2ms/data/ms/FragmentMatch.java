package org.ms2ms.data.ms;

public class FragmentMatch implements Comparable<FragmentMatch>
{
  private FragmentEntry mEntry;
  private int           mObsIdx;

  public FragmentMatch() { super(); }
  public FragmentMatch(FragmentEntry s, int i)
  {
    mEntry=s; mObsIdx=i;
  }

  public FragmentEntry getEntry()    { return mEntry; }
  public int           getObsIndex() { return mObsIdx; }

  @Override
  public int compareTo(FragmentMatch o) { return mEntry!=null?mEntry.compareTo(o.mEntry):0; }
}
