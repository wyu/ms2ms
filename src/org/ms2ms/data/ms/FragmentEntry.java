package org.ms2ms.data.ms;

import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

/**
 * Created by yuw on 8/4/16.
 */
public class FragmentEntry implements Comparable<FragmentEntry>
{
  private int   mLength;
  private long  mPeptideKey;
  private float mMH;
  private FragmentEntry mPrev;

  public FragmentEntry() { super(); }
  public FragmentEntry(FragmentEntry s)
  {
    super();
    mLength=s.mLength; mPeptideKey=s.mPeptideKey; mMH=s.mMH; mPrev=s.mPrev;
  }
  public FragmentEntry(Float mh, Long peptide, FragmentEntry prev, int len)
  {
    super();
    mMH=mh; mPeptideKey=peptide; mPrev=prev; mLength=len;
  }

  public int           getLength()     { return mLength; }
  public Float         getMH()         { return mMH; }
  public Long          getPeptideKey() { return mPeptideKey; }
  public FragmentEntry getPrev()       { return mPrev; }

  public FragmentEntry setLen(int s)   { mLength=s; return this; }

  public FragmentEntry copy()
  {
    return new FragmentEntry(mMH, mPeptideKey, mPrev, mLength);
  }
  @Override
  public int compareTo(FragmentEntry o)
  {
    int c = Long.compare(mPeptideKey, o.getPeptideKey());
    if (c==0) c = Integer.compare(mLength, o.getLength());

    return c;
  }

  @Override
  public String toString()
  {
    return mPeptideKey+":"+mLength+"#"+ (Tools.d2s(mMH, 4))+"<-"+(mPrev!=null?Tools.d2s(mPrev.getMH(), 4):"NUL");
  }
}
