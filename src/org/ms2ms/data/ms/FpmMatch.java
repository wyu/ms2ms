package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.ms2ms.utils.Tools;

public class FpmMatch implements Comparable<FpmMatch>
{
  private FpmEntry mEntry;
  private int      mProteinIdx, mPeptideSeqIdx;
  private IonType  mIonType;

  public FpmMatch() { super(); }
  public FpmMatch(FpmEntry s, IonType ion, Integer protein, Integer seq)
  {
    mEntry=s; mProteinIdx=protein; mIonType=ion; mPeptideSeqIdx=seq;
  }

  public FpmEntry getEntry()      { return mEntry; }
  public Integer  getProtein()    { return mProteinIdx; }
  public IonType  getIonType()    { return mIonType; }
  public int      getPeptideSeq() { return mPeptideSeqIdx; }

  // order by the peptide key in desending order
  @Override
  public int compareTo(FpmMatch o)
  {
    int c = Integer.compare(o.mProteinIdx, mProteinIdx);
//    if (c==0) c = Integer.compare(
//        o.getEntry().getTrack().size(),
//          getEntry().getTrack().size());
    if (c==0) c =   Double.compare(o.getEntry().getGapScore(), getEntry().getGapScore());
    if (c==0) c =  Integer.compare(o.mPeptideSeqIdx, mPeptideSeqIdx);
    if (c==0) c = mIonType.compareTo(o.mIonType);

    return c;
  }

  @Override
  public String toString()
  {
    return mIonType.toString()+",protein:"+mProteinIdx+"; "+mEntry.toString();
  }

}

