package org.ms2ms.data.ms;

import org.ms2ms.Disposable;
import org.ms2ms.data.Binary;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Tools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

/**
 * Slim-down version of FragmentEntry, Created by yuw on 7/8/2019.
 */
public class FragmentSlot implements Comparable<FragmentSlot>, Disposable, Binary
{
  protected byte  mCharge=1;
  protected byte  mLength;
  protected int   mPeptideKey;
  protected float mMH;

  public FragmentSlot() { super(); }
  public FragmentSlot(FragmentEntry s)
  {
    super();
    if (s!=null) { mLength=(byte )s.getLength(); mPeptideKey=s.getPeptideKey(); mMH=s.getMH(); }
  }
  public FragmentSlot(FragmentSlot s)
  {
    super();
    if (s!=null) { mLength=s.mLength; mPeptideKey=s.mPeptideKey; mMH=s.mMH; }
  }
  @Deprecated
  public FragmentSlot(Float mh, int peptide, int len)
  {
    super();
    mMH=mh; mPeptideKey=peptide; mLength=(byte )len;
  }
  public FragmentSlot(float mh, int z, int peptide, int len)
  {
    super();
    mMH=mh; mPeptideKey=peptide; mLength=(byte )len; mCharge=(byte )z;
  }

  public int           getLength()       { return Math.abs(mLength); }
  public int           getCharge()       { return mCharge; }
  public float         getMH()           { return mMH; }
  public int           getPeptideKey()   { return mPeptideKey; }

  public FragmentSlot isProDirected(boolean s)
  {
    if (s) mLength = (byte )(Math.abs(mLength)*-1);
    return this;
  }
  public boolean isProDirected() { return mLength<0; }
  public FragmentSlot setCharge(int s) { mCharge=(byte )s; return this; }
  public FragmentSlot setPeptideKey(Integer s)
  {
    mPeptideKey=s;
    return this;
  }
  public FragmentSlot setLen(int s)   { mLength=(byte )s; return this; }
  public FragmentSlot setMH(float s)   { mMH    =s; return this; }

  public FragmentSlot copy()
  {
    return new FragmentSlot(mMH, mCharge, mPeptideKey, mLength);
  }
  @Override
  public int compareTo(FragmentSlot o)
  {
    int c = Long.compare(mPeptideKey, o.getPeptideKey());
    if (c==0) c = Integer.compare(getLength(), o.getLength());

    return c;
  }

  public int compareByPeptide(FragmentSlot o)
  {
    return Long.compare(mPeptideKey, o.getPeptideKey());
  }

  @Override
  public String toString()
  {
    return mPeptideKey+":"+getLength()+"#"+ (Tools.d2s(mMH, 4));
  }

  @Override
  public void dispose() { }

  @Override
  public void write(DataOutput ds) throws IOException
  {
    IOs.write(ds, mLength);
    IOs.write(ds, mCharge);
    IOs.write(ds, mPeptideKey);
    IOs.write(ds, mMH);
  }

  @Override
  public void read(DataInput ds) throws IOException
  {
    mLength = IOs.read(ds, mLength);
    mCharge = IOs.read(ds, mCharge);
    mPeptideKey = IOs.read(ds, mPeptideKey);
    mMH = IOs.read(ds, mMH);
  }
}
