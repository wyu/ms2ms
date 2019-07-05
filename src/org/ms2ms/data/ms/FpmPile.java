package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.spectrum.IonType;

public class FpmPile extends AbstractPile<FpmMatch>
{
  private ProteinSegments mProtein = new ProteinSegments(IonType.y, IonType.b);

  public FpmPile()       { super(); }
  public FpmPile(int s)  { super(s); init(); }

  @Override
  public int getKeyAt(int pile, int idx) { return get(pile,idx).getProtein(); }
  @Override
  public FpmMatch[] newPile(int s) { return new FpmMatch[s]; }

  public ProteinSegments getProtein() { return mProtein; }
  public ProteinSegments nextProtein()
  {
    // once a block of FragmentEntry was located in the 'pile', we can turn them into a 'track'
    mProtein.ensureCapacity(mSeriesEnd, 1.5).start();

    // skip the repeats
    double scr=-1000; int seq=0;
    for (int i=0; i<mSeriesEnd; i++)
    {
      if (at(i).getPeptideSeq()!=seq || at(i).getEntry().getGapScore()!=scr)
      {
        mProtein.addSegment(at(i));
        scr = at(i).getEntry().getGapScore(); seq = at(i).getPeptideSeq();
      }
    }

    return mProtein;
  }
}