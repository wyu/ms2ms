package org.ms2ms.data.ms;

import java.util.Collection;

/**
 * Created by yuw on 3/31/17.
 */
public class FragmentEntries
{
  private FragmentEntry[] mEntries;

  public FragmentEntries() { super(); }
  public FragmentEntries(Collection<FragmentEntry> s)
  {
    if (s!=null && s.size()>0)
    {
      mEntries = new FragmentEntry[s.size()]; int idx=0;
      for (FragmentEntry E : s) mEntries[idx++]=E;
    }
    mEntries=null;
  }
}
