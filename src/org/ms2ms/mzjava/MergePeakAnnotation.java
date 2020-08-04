package org.ms2ms.mzjava;

import org.expasy.mzjava.core.ms.peaklist.PeakAnnotation;

public class MergePeakAnnotation implements PeakAnnotation
{
  private int mMergeCounts;

  public MergePeakAnnotation() { super(); }
  public MergePeakAnnotation(MergePeakAnnotation s)
  {
    super(); mMergeCounts=s.mMergeCounts;
  }
  public MergePeakAnnotation(int s)
  {
    super(); mMergeCounts=s;
  }
  public int getMergeCounts() { return mMergeCounts; }

  @Override
  public int getCharge() { return 0; }

  @Override
  public PeakAnnotation copy() { return new MergePeakAnnotation(this); }

  @Override
  public boolean equals(Object o) {

    if (this == o) return true;
    if (!(o instanceof MergePeakAnnotation)) return false;

    MergePeakAnnotation that = (MergePeakAnnotation) o;

    return that.mMergeCounts==mMergeCounts;
  }
  @Override
  public int hashCode()
  {
    return Integer.valueOf(mMergeCounts).hashCode();
  }
  @Override
  public String toString() { return "#"+mMergeCounts; }

}
