package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.spectrum.IonType;

import java.util.Comparator;

public class FpmMatch implements Comparable<FpmMatch>
{
  private FpmSlot mEntry;
  private int      mProteinIdx, mPeptideSeqIdx;
  private IonType  mIonType;

  public static class TrackDesendComparator implements Comparator<FpmMatch>
  {
    public TrackDesendComparator() { }
    public int compare(FpmMatch o1, FpmMatch o2)
    {
      if (o2.getEntry()==null) return -1;
      if (o1.getEntry()==null) return  1;

      int c = Integer.compare(o2.getEntry().getTrack().size(), o1.getEntry().getTrack().size());
      if (c==0) c = Integer.compare(o2.getEntry().getMotifs(), o1.getEntry().getMotifs());

      return c;
    }
  }
  public static class ProteinDesendComparator implements Comparator<FpmMatch>
  {
    public ProteinDesendComparator() { }
    public int compare(FpmMatch o1, FpmMatch o2)
    {
      if (o2.getEntry()==null) return -1;
      if (o1.getEntry()==null) return  1;

      int c = Float.compare(Math.signum(o2.getProtein()), Math.signum(o1.getProtein()));
      if (c==0) c = Integer.compare(o1.getProtein(), o2.getProtein());
      if (c==0) c = Integer.compare(o2.getEntry().getTrack().size(), o1.getEntry().getTrack().size());

      return c;
    }
  }

  public FpmMatch() { super(); }
  public FpmMatch(FpmSlot s, IonType ion, int protein, int seq)
  {
    mEntry=s; mProteinIdx=protein; mIonType=ion; mPeptideSeqIdx=seq;
  }

  public FpmSlot getEntry()       { return mEntry; }
  public int      getProtein()    { return mProteinIdx; }
  public IonType  getIonType()    { return mIonType; }
  public int      getPeptideSeq() { return mPeptideSeqIdx; }

  // order by the peptide key in desending order
  @Override
  public int compareTo(FpmMatch o)
  {
    int c = Integer.compare(o.mProteinIdx, mProteinIdx);
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

