package org.ms2ms.data.ms;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.UnmodifiableIterator;
import com.sun.xml.bind.v2.schemagen.xmlschema.Annotated;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.ms2ms.Disposable;
import org.ms2ms.algo.Peaks;
import org.ms2ms.math.Stats;
import org.ms2ms.mzjava.AnnotatedPeak;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.util.*;

/** The Fragment(predicted)-peak(obs)-match entry
 *
 * Created by yuw on 8/6/16.
 */
public class FpmEntry implements Comparable<FpmEntry>, Disposable
{
  private boolean mHas1st=false, mExpectedY1=false;
  private int                 mMotifs=0, m1stPass=0, mWeaks=0, mC13=0;
  private double              mGapScore=0, mIntensities=0d;
  private Double              mMatchScore=null;
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
      mGapScore=f.mGapScore; mIntensities=f.mIntensities;
//      mProb=f.mProb; mKaiScore=f.mKaiScore;
      mFragment = new FragmentEntry(f.mFragment);
      mTrack = f.mTrack.asList();
    }
  }
  public FpmEntry(FragmentEntry f, List<AnnotatedPeak> t)
  {
    super();
    mFragment=f; mTrack=ImmutableList.copyOf(t);
  }
//  public FpmEntry(FragmentEntry f, List<AnnotatedPeak> t, int motifs, int size_1st)
//  {
//    super();
//    mFragment=f; mTrack=ImmutableList.copyOf(t); mMotifs=motifs; m1stPass=size_1st;
//  }
//
//  public FpmEntry(FragmentEntry f, List<AnnotatedPeak> t, double ai)
//  {
//    super();
//    mFragment=f; mTrack=ImmutableList.copyOf(t); mIntensities=ai;
//  }

  public boolean             isDecoy()      { return mFragment!=null && mFragment.getPeptideKey()!=null && mFragment.getPeptideKey()<0; }
  public boolean             has1st()       { return mHas1st; }
//  public double              getKaiScore()  { return mKaiScore; }
  public double              getGapScore()  { return mGapScore; }
  public double              getIntensity() { return mIntensities; }
//  public double              getProb()      { return mProb; }
  public int                 getMotifs()    { return mMotifs; }
  public Double             getMatchScore() { return mMatchScore; }
//  public int                 get1stPass()   { return m1stPass; }
  public FragmentEntry       getFragment()  { return mFragment; }
  public ImmutableList<AnnotatedPeak> getTrack()     { return mTrack; }
  public AnnotatedPeak       at(int s)      { return mTrack.get(s); }

//  public FpmEntry increIntensities(double s) { mIntensities+=s; return this; }

  public FpmEntry has1st(        boolean s) { mHas1st    =s;  return this; }
  public FpmEntry hasExpectedY1( boolean s) { mExpectedY1=s;  return this; }
  public FpmEntry setMotifs(         int s) { mMotifs    =s;  return this; }
  public FpmEntry setMatchScore(  Double s) { mMatchScore=s; return this; }

  //  public FpmEntry set1stPass(        int s) { m1stPass   =s; return this; }
//  public FpmEntry setWeaks(          int s) { mWeaks     =s; return this; }
//  public FpmEntry setC13(            int s) { mC13       =s; return this; }
//  public FpmEntry setProb(        double s) { mProb      =s; return this; }
  public FpmEntry setGapScore(    double s) { mGapScore  =s; return this; }
//  public FpmEntry setStdErrRegression(double s) { mKaiScore  =s; return this; }
  public FpmEntry setIntensity(   double s) { mIntensities=s; return this; }

  public FpmEntry shallow_copy()
  {
    FpmEntry clone = new FpmEntry();
    // no scores, just the skeleton
    clone.mFragment = mFragment;
    clone.mTrack=mTrack;

    return clone;
  }

  // update the ion matches and return the local base intensity
  public double match(OffsetPpmTolerance tol)
  {
    double local_base=0;
    if (Tools.isSet(getTrack()))
    {
      List<AnnotatedPeak> pts = new ArrayList<>(getTrack().size());
      for (AnnotatedPeak  pt : getTrack())
      {
        Double calc = pt.getProperty("calc"), obs = pt.getIntensity();
        if (calc!=null && obs!=null && tol.withinTolerance(calc, obs))
        {
          AnnotatedPeak p = pt.clone();
          p.setMzAndCharge(Stats.ppm(obs, calc)+tol.getOffset(obs), pt.getCharge());
          pts.add(p);
          if (p.getSNR()>local_base) local_base=p.getSNR();
        }
      }
      mTrack = ImmutableList.copyOf(pts);
    }

    return local_base;
  }
  // following the re-match, update the scores of the series
  public FpmEntry score(Double local_base)
  {
    if (Tools.isSet(getTrack()))
    {
      // we used to score only the y1.
      int best=0, start=1, delta=0;
      double percentile=0, score=(has1st()?calcGapScore(getTrack().get(getTrack().size()-1), 1, 1d)*0.1d:0),
          sumAI=0d, base=local_base!=null?1d/local_base:0.01d;
      for (int i=getTrack().size()-1; i>=0; i--)
      {
        AnnotatedPeak pk = getTrack().get(i);

        delta = pk.getCharge()-start;
        percentile = (pk.getSNR()*base); // set a minimum
        // accumualte the gap score
        if (delta>0) score+=(calcGapScore(pk, delta, 1d))*percentile;

        start=pk.getCharge(); sumAI+=percentile;

        int first=getTrack().get(i).getCharge(), last=first;
        for (int j=i+1; j<getTrack().size(); j++)
          if (last-getTrack().get(j).getCharge()!=1) break; else last=getTrack().get(j).getCharge();

        if (first-last>best) best=first-last;
      }
      setMotifs(best).setGapScore(-10d*score);
    }

    return this;
  }
  // BEWARE that the information are saved in the 'match' object in a non-std way.
  private double calcGapScore(AnnotatedPeak match, int gap, double min_ppm)
  {
    // no point to proceed...
    if (gap<=0/* || gap>4*/) return 0;

    // average AA mass: 115, largest - smallest: 129
    int     bins = (int )(Math.log(2)/Math.log(1d+1E-6*Math.max(Math.abs(match.getMz()), min_ppm))),
        nsamples = 0, ntrials=(int )Math.round(((gap-1)*115d+129d)*match.getFrequency());

    // cumulative numbers of gaps
    for (int i=1; i<=gap; i++)
      if (i<19 && nsamples<bins) nsamples+=Math.exp(Stats.ln_combination(19, i)); else break;

    double score0=0;
    if (nsamples<bins/2)
    {
      score0 = -0.07491232 + 0.41163668*Math.log((double )nsamples/(double )bins) + 0.40504996*Math.log((double )ntrials);
      if (score0 > -0.1) score0=0;
    }
    match.setProperty("score", score0);

    return score0;
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
    return "#"+(Tools.isSet(mTrack)?mMotifs+"/"+mTrack.size():0)+"@"+(mFragment!=null?mFragment.toString():"");
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
  public int hashcodeByTrackMz(int precision)
  {
    int hash=0; double multi = Math.pow(10d, precision);
    if (Tools.isSet(getTrack()))
      for (int i=0; i<getTrack().size(); i++)
        // grab the obs m/z rather than the ppm
        hash+=Math.round(getTrack().get(i).getIntensity()*multi);

    return hash;
  }
  public Collection<Double> fillMz(Collection<Double> data)
  {
    if (Tools.isSet(getTrack()))
      for (AnnotatedPeak p : getTrack()) data.add(p.getIntensity());

    return data;
  }
  public FpmEntry inspect(Set<Double> y1s, Double local_base)
  {
    if (!Tools.isSet(getTrack())) return this;

    boolean has1st = false, y1=false;
    AnnotatedPeak pk1 = null;
    List<AnnotatedPeak> track = new ArrayList<>(getTrack().size());
    for (int i=0; i<getTrack().size(); i++)
    {
      track.add(at(i));
      // check the presence of y1
      if (at(i).getCharge()==1 && !has1st) { has1st=true; pk1=at(i); }
      // the observed m/z was saved in the 'intensity' field!!!
      if (at(i).getCharge()==1 && y1s!=null && y1s.contains(at(i).getIntensity())) y1=true;
    }

    // return now if there isn't enough matches
    if (track.size()==0)
      return has1st(has1st).hasExpectedY1(Tools.isSet(y1s) && y1).setIntensity(0).setMotifs(0).setGapScore(0);

    // we used to score only the y1.
    int best=0, start=1, delta=0;
    double percentile=0, score=(has1st?calcGapScore(track.get(track.size() - 1), 1, 1d)*0.1d:0),
        sumAI=0d, base=local_base!=null?1d/local_base:0.01d;
    for (int i=track.size()-1; i>=0; i--)
    {
      AnnotatedPeak pk = track.get(i);

      delta = pk.getCharge()-start;
      percentile = (pk.getSNR()*base); // set a minimum
      // accumualte the gap score
      if (delta>0) score+=(calcGapScore(pk, delta, 1d))*percentile;

      start=pk.getCharge(); sumAI+=percentile;

      int first=track.get(i).getCharge(), last=first;
      for (int j=i+1; j<track.size(); j++)
        if (last-track.get(j).getCharge()!=1) break; else last=track.get(j).getCharge();

      if (first-last>best) best=first-last;
    }

    // http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/tutorials/xegbohtmlnode16.html
    return has1st(has1st).hasExpectedY1(Tools.isSet(y1s) && y1).setIntensity(sumAI).setMotifs(best).setGapScore(-10d*score);
  }

  @Override
  public FpmEntry clone()
  {
    return new FpmEntry(this);
  }
  @Override
  public int hashCode()
  {
    int hash = Tools.hashCodes(mHas1st,mExpectedY1);
        hash+= mMotifs+m1stPass+mWeaks+mC13;
        hash+= Tools.hashCodes(mIntensities,mGapScore);
//        hash+= Tools.hashCodes(mIntensities,mProb,mGapScore,mKaiScore);
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
