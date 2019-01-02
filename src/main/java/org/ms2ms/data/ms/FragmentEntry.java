package org.ms2ms.data.ms;

import org.ms2ms.Disposable;
import org.ms2ms.data.Binary;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Tools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

/**
 * Created by yuw on 8/4/16.
 */
public class FragmentEntry implements Comparable<FragmentEntry>, Disposable, Binary
{
  private char  mCharge=1;
  private int   mLength;
  private int   mPeptideKey;
  private float mMH;
  private FragmentEntry mPrev;

  public FragmentEntry() { super(); }
  public FragmentEntry(FragmentEntry s)
  {
    super();
    mLength=s.mLength; mPeptideKey=s.mPeptideKey; mMH=s.mMH; mPrev=s.mPrev;
  }
  @Deprecated
  public FragmentEntry(Float mh, int peptide, FragmentEntry prev, int len)
  {
    super();
    mMH=mh; mPeptideKey=peptide; mPrev=prev; mLength=len;
  }
  public FragmentEntry(Float mh, char z, int peptide, FragmentEntry prev, int len)
  {
    super();
    mMH=mh; mPeptideKey=peptide; mPrev=prev; mLength=len; mCharge=z;
  }

  public int           getLength()       { return mLength; }
  public char          getCharge()       { return mCharge; }
  public float         getMH()           { return mMH; }
  public int           getPeptideKey()   { return mPeptideKey; }
  public FragmentEntry getPrev()         { return mPrev; }

  public FragmentEntry setCharge(char s) { mCharge=s; return this; }
  public FragmentEntry setPeptideKey(Integer s)
  {
    mPeptideKey=s;
    return this;
  }
  public FragmentEntry setLen(           int s)   { mLength=s; return this; }
  public FragmentEntry setMH(          float s)   { mMH    =s; return this; }
  public FragmentEntry setPrev(FragmentEntry s)   { mPrev  =s; return this; }

  public FragmentEntry copy()
  {
    return new FragmentEntry(mMH, mCharge, mPeptideKey, mPrev, mLength);
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

  @Override
  public void dispose() { mPrev=null; }

  @Override
  public void write(DataOutput ds) throws IOException
  {
    IOs.write(ds, mLength);
    IOs.write(ds, mCharge);
    IOs.write(ds, mPeptideKey);
    IOs.write(ds, mMH);
    // can't write the actual object. to avoid recursive
    // write the peptide seq key and MH so we can hook up the right frag at later time
    IOs.write(ds, mPrev!=null?mPrev.mMH:0f);
  }

  @Override
  public void read(DataInput ds) throws IOException
  {
    mLength = IOs.read(ds, mLength);
    mCharge = IOs.read(ds, mCharge);
    mPeptideKey = IOs.read(ds, mPeptideKey);
    mMH = IOs.read(ds, mMH);

    float p = IOs.read(ds, mMH);
    if (p!=0f) mPrev = new FragmentEntry(p, '0', 0, null, 0);
  }
}
