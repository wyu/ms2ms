package org.ms2ms.data.ms;

public class FragmentMatch implements Comparable<FragmentMatch>
{
  private FragmentSlot mEntry;
  private short        mObsIdx; // -32,767 to 32,767, 2 bytes

  public FragmentMatch() { super(); }
  public FragmentMatch(FragmentSlot s, int i)
  {
    mEntry=s; mObsIdx=(short )i;
  }

  public FragmentSlot getEntry()    { return mEntry; }
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
