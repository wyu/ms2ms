package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.peaklist.Peak;

import java.util.List;

/** The Fragment(predicted)-peak(obs)-match entry
 *
 * Created by yuw on 8/6/16.
 */
public class FpmEntry implements Comparable<FpmEntry>
{
  private Double        mIntensities=0d;
  private FragmentEntry mFragment   =null;
  private List<Peak>    mTrack      =null;

  public FpmEntry() { super(); }
  public FpmEntry(FragmentEntry f, List<Peak> t)
  {
    super();
    mFragment=f; mTrack=t;
  }
  public FpmEntry(FragmentEntry f, List<Peak> t, double ai)
  {
    super();
    mFragment=f; mTrack=t; mIntensities=ai;
  }

  public FragmentEntry getFragment() { return mFragment; }
  public List<Peak>    getTrack()    { return mTrack; }

  public FpmEntry increIntensities(double s) { mIntensities+=s; return this; }

  @Override
  public int compareTo(FpmEntry o)
  {
    int c = mFragment.compareTo(o.getFragment());
    if (c==0 && mTrack!=null && o.getTrack()!=null) c = Integer.compare(mTrack.size(), o.getTrack().size());

    return c;
  }
}
