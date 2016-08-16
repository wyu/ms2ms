package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.ms2ms.mzjava.AnnotatedPeak;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/** The Fragment(predicted)-peak(obs)-match entry
 *
 * Created by yuw on 8/6/16.
 */
public class FpmEntry implements Comparable<FpmEntry>
{
  private boolean mHas1st=false, mExpectedY1=false;
  private int                 mMotifs=0, m1stPass=0, mWeaks=0, mC13=0;
  private double              mIntensities=0d, mProb=0, mGapScore=0, mKaiScore=0;
  private FragmentEntry       mFragment   =null;
  private List<AnnotatedPeak> mTrack      =null;

  public FpmEntry() { super(); }
  public FpmEntry(FragmentEntry f)
  {
    super();
    mFragment=f;
  }
  public FpmEntry(FragmentEntry f, List<AnnotatedPeak> t)
  {
    super();
    mFragment=f; mTrack=t;
  }
  public FpmEntry(FragmentEntry f, List<AnnotatedPeak> t, int motifs, int size_1st)
  {
    super();
    mFragment=f; mTrack=t; mMotifs=motifs; m1stPass=size_1st;
  }

  public FpmEntry(FragmentEntry f, List<AnnotatedPeak> t, double ai)
  {
    super();
    mFragment=f; mTrack=t; mIntensities=ai;
  }

  public double              getKaiScore()  { return mKaiScore; }
  public double              getGapScore()  { return mGapScore; }
  public double              getIntensity() { return mIntensities; }
  public double              getProb()      { return mProb; }
  public int                 getMotifs()    { return mMotifs; }
  public int                 get1stPass()   { return m1stPass; }
  public FragmentEntry       getFragment()  { return mFragment; }
  public List<AnnotatedPeak> getTrack()     { return mTrack; }
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
  public FpmEntry setKaiScore(    double s) { mKaiScore  =s; return this; }
  public FpmEntry setIntensity(   double s) { mIntensities=s; return this; }

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
    return "#"+(Tools.isSet(mTrack)?mMotifs+"/"+mTrack.size()+"/"+m1stPass:0)+"@"+mFragment.toString();
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
}
