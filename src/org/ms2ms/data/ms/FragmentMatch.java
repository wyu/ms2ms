package org.ms2ms.data.ms;

import org.ms2ms.utils.Tools;

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

  // order by the peptide key in desending order
  @Override
  public int compareTo(FragmentMatch o) { return mEntry!=null?o.mEntry.compareTo(mEntry):0; }

  @Override
  public String toString()
  {
    return "#"+mObsIdx+"; "+mEntry.toString();
  }

}
