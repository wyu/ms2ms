package org.ms2ms.data.ms;

import com.google.common.collect.ImmutableList;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.ms2ms.Disposable;
import org.ms2ms.mzjava.AnnotatedPeak;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/** The Fragment(predicted)-peak(obs)-match entry
 *
 * Created by yuw on 8/6/16.
 */
public class FpmEntry implements Comparable<FpmEntry>, Disposable
{
  private boolean mHas1st=false, mExpectedY1=false;
  private int                 mMotifs=0, m1stPass=0, mWeaks=0, mC13=0;
  private double              mIntensities=0d, mProb=0, mGapScore=0, mKaiScore=0;
  private FragmentEntry       mFragment   =null;
  private ImmutableList<AnnotatedPeak> mTrack      =null;

  public FpmEntry()
  {
    super(); mFragment = new FragmentEntry();
  }
  public FpmEntry(FpmEntry f)
  {
    super();
    if (f!=null)
    {
      mHas1st=f.mHas1st; mExpectedY1=f.mExpectedY1; mMotifs=f.mMotifs; m1stPass=f.m1stPass; mWeaks=f.mWeaks; mC13=f.mC13;
      mIntensities=f.mIntensities; mProb=f.mProb; mGapScore=f.mGapScore; mKaiScore=f.mKaiScore;
      mFragment = new FragmentEntry(f.mFragment);
      mTrack=f.mTrack;
    }
  }
  public FpmEntry(FragmentEntry f)
  {
    super();
    mFragment=f;
  }
  public FpmEntry(FragmentEntry f, List<AnnotatedPeak> t)
  {
    super();
    mFragment=f; mTrack=ImmutableList.copyOf(t);
  }
  public FpmEntry(FragmentEntry f, List<AnnotatedPeak> t, int motifs, int size_1st)
  {
    super();
    mFragment=f; mTrack=ImmutableList.copyOf(t); mMotifs=motifs; m1stPass=size_1st;
  }

  public FpmEntry(FragmentEntry f, List<AnnotatedPeak> t, double ai)
  {
    super();
    mFragment=f; mTrack=ImmutableList.copyOf(t); mIntensities=ai;
  }

  public boolean             isDecoy()      { return mFragment!=null && mFragment.getPeptideKey()!=null && mFragment.getPeptideKey()<0; }
  public boolean             has1st()       { return mHas1st; }
  public double              getKaiScore()  { return mKaiScore; }
  public double              getGapScore()  { return mGapScore; }
  public double              getIntensity() { return mIntensities; }
  public double              getProb()      { return mProb; }
  public int                 getMotifs()    { return mMotifs; }
  public int                 get1stPass()   { return m1stPass; }
  public FragmentEntry       getFragment()  { return mFragment; }
  public ImmutableList<AnnotatedPeak> getTrack()     { return mTrack; }
  public AnnotatedPeak       at(int s)      { return mTrack.get(s); }

  public FpmEntry increIntensities(double s) { mIntensities+=s; return this; }

  public FpmEntry has1st(        boolean s) { mHas1st    =s; return this; }
  public FpmEntry hasExpectedY1( boolean s) { mExpectedY1=s; return this; }
  public FpmEntry setMotifs(         int s) { mMotifs    =s; return this; }
  public FpmEntry set1stPass(        int s) { m1stPass   =s; return this; }
  public FpmEntry setWeaks(          int s) { mWeaks     =s; return this; }
  public FpmEntry setC13(            int s) { mC13       =s; return this; }
  public FpmEntry setProb(        double s) { mProb      =s; return this; }
  public FpmEntry setGapScore(    double s) { mGapScore  =s; return this; }
  public FpmEntry setStdErrRegression(double s) { mKaiScore  =s; return this; }
  public FpmEntry setIntensity(   double s) { mIntensities=s; return this; }

  public FpmEntry shallow_copy()
  {
    FpmEntry clone = new FpmEntry();
    // no scores, just the skeleton
    clone.mFragment = mFragment;
    clone.mTrack=mTrack;

    return clone;
  }

  @Override
  public int compareTo(FpmEntry o)
  {
    int c = mFragment.compareTo(o.getFragment());
    if (c==0 && mTrack!=null && o.getTrack()!=null) c = Integer.compare(mTrack.size(), o.getTrack().size());

    return c;
  }

  @Override
  public String toString()
  {
    return "#"+(Tools.isSet(mTrack)?mMotifs+"/"+mTrack.size():0)+"@"+mFragment.toString();
  }
  public static Map<String, Object> report(List<AnnotatedPeak> track)
  {
    Map<String, Object> stats = null;

    String ppms=null, pos=null; double ais=0;
    if (Tools.isSet(track))
    {
      for (Peak pk : track)
      {
        ppms = Strs.extend(ppms, Tools.d2s(pk.getMz(), 1), ",");
        pos  = Strs.extend(pos,  pk.getCharge()+"",        ",");
        ais += Math.log(pk.getIntensity());
      }
      stats = new HashMap<>();
      stats.put("ppm track", ppms);
      stats.put("position track", pos);
      stats.put("ai logSum", ais);
    }

    return stats;
  }
  public int hashcodeByTrack()
  {
    int hash=0;
    if (Tools.isSet(getTrack()))
      for (int i=0; i<getTrack().size(); i++)
        hash+=(i+1)*getTrack().get(i).getMz()*1000;

    return hash;
  }
  public int hashcodeByTrackMz()
  {
    int hash=0;
    if (Tools.isSet(getTrack()))
      for (int i=0; i<getTrack().size(); i++)
        hash+=getTrack().get(i).getMz()*1000;

    return hash;
  }
  @Override
  public int hashCode()
  {
    int hash = Tools.hashCodes(mHas1st,mExpectedY1);
        hash+= mMotifs+m1stPass+mWeaks+mC13;
        hash+= Tools.hashCodes(mIntensities,mProb,mGapScore,mKaiScore);
    if (mFragment==null) hash+=mFragment.hashCode();

    return hash+hashcodeByTrack();
  }

  @Override
  public void dispose()
  {
    Tools.dispose(mFragment);
    mTrack=null;
  }
}
