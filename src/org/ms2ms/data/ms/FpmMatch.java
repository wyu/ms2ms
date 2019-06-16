package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.spectrum.IonType;

public class FpmMatch implements Comparable<FpmMatch>
{
  private FpmEntry mEntry;
  private int      mProteinIdx;
  private IonType  mIonType;

  public FpmMatch() { super(); }
  public FpmMatch(FpmEntry s, IonType ion, Integer protein)
  {
    mEntry=s; mProteinIdx=protein; mIonType=ion;
  }

  public FpmEntry getEntry()   { return mEntry; }
  public Integer  getProtein() { return mProteinIdx; }
  public IonType  getIonType() { return mIonType; }

  // order by the peptide key in desending order
  @Override
  public int compareTo(FpmMatch o)
  {
    int c = Integer.compare(o.mProteinIdx, mProteinIdx);
    if (c==0) c = Integer.compare(
        o.getEntry().getTrack().size(),
          getEntry().getTrack().size());

    return c;
  }
}

