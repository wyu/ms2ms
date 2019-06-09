package org.ms2ms.data.ms;

import com.google.common.collect.ImmutableList;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.ms2ms.Disposable;
import org.ms2ms.data.Binary;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.*;

/** The Fragment(predicted)-peak(obs)-match entry
 *
 * Created by yuw on 8/6/16.
 */
public class FpmEntry implements Comparable<FpmEntry>, Disposable, Binary
{
  private boolean mHas1st=false, mExpectedY1=false;
  private int                 mMotifs=0, m1stPass=0, mWeaks=0, mC13=0, mPros=0;
  private double              mGapScore=0, mIntensities=0d, mGapScore0=0;
//  private Double              mMatchScore =null;
  private FragmentEntry       mFragment   =null;
  private ImmutableList<PeakMatch> mTrack =null;

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
      mFragment = new FragmentEntry(f.mFragment);
      mTrack = f.mTrack!=null?f.mTrack.asList():null;
    }
  }
  public FpmEntry(FragmentEntry f, List<PeakMatch> t)
  {
    super();
    mFragment=f; mTrack=ImmutableList.copyOf(t);
  }
  public FpmEntry(FragmentEntry f, PeakMatch[] matches, int size)
  {
    super();
    mFragment=f; mTrack=ImmutableList.copyOf(Arrays.copyOfRange(matches, 0, size));
  }

  public boolean             isDecoy()      { return mFragment!=null && mFragment.getPeptideKey()<0; }
  public boolean             has1st()       { return mHas1st; }
  public boolean             hasY1()        { return mExpectedY1; }
  public double              getGapScore()  { return mGapScore; }
  public double              getIntensity() { return mIntensities; }
  public int                 getMotifs()    { return mMotifs; }
  public int                 getProlines()  { return mPros; }
  public FpmEntry            increProlines()  { mPros++; return this; }
//  public Double             getZScore() { return mMatchScore; }
  public FragmentEntry       getFragment()  { return mFragment; }
  public ImmutableList<PeakMatch> getTrack()     { return mTrack; }
  public PeakMatch       at(int s)      { return mTrack.get(s); }

  public FpmEntry has1st(        boolean s) { mHas1st    =s;  return this; }
  public FpmEntry hasExpectedY1( boolean s) { mExpectedY1=s;  return this; }
  public FpmEntry setMotifs(         int s) { mMotifs    =s;  return this; }

  public FpmEntry setGapScore(    double s) { mGapScore  =s; return this; }
  public FpmEntry setGapScore0(   double s) { mGapScore0 =s; return this; }
  public FpmEntry setIntensity(   double s) { mIntensities=s; return this; }
  public FpmEntry setTrack(ImmutableList<PeakMatch> s) { mTrack=s; return this; }

  // need the calc mz of the unmodified backbone in the order of the ion series
  public Long[] hashByWTCalcMz(double[] backbone)
  {
    if (!Tools.isSet(getTrack())) return null;

    // by the calculated mass and position
    int[] ions = new int[backbone.length]; int stop=-1;
    for (int i=0; i<getTrack().size(); i++)
    {
      int k=getTrack().get(i).getCharge()-1;
      if (k>=backbone.length)
        continue;
      ions[k] = getTrack().size()-1-i;
      if (k>stop) stop=k;
    }

    Long   hash=0L; double m=0d; int k=0;
    Long[] hashes = new Long[getTrack().size()];
    for (int i=0; i<=stop; i++)
    {
      m+=backbone[i]; hash += (i+1)*Double.hashCode(m);
      if (ions[i]>0) hashes[ions[i]]=hash;
    }

    return hashes;
  }
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
      List<PeakMatch> pts = new ArrayList<>(getTrack().size());
      for (PeakMatch  pt : getTrack())
      {
        Double calc = pt.getScore(), obs = pt.getIntensity();
        if (calc!=null && obs!=null && tol.withinTolerance(calc, obs))
        {
          PeakMatch p = pt.clone();
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
        PeakMatch pk = getTrack().get(i);

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
  private double calcGapScore(PeakMatch match, int gap, double min_ppm)
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
    match.setScore(score0);

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
    return "#"+(Tools.isSet(mTrack)?mMotifs+"/"+mTrack.size():0)+"@"+(mFragment!=null?mFragment.toString():"")+"$"+Tools.d2s(getGapScore(), 1);
  }
  public static Map<String, Object> report(List<PeakMatch> track)
  {
    Map<String, Object> stats = null;

    String ppms=null, pos=null; double ais=0;
    if (Tools.isSet(track))
    {
      for (PeakMatch pk : track)
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
      for (PeakMatch p : getTrack()) data.add(p.getIntensity());

    return data;
  }
  // estimates the preliminary score for the screening step only!
  // NO re-calibration or intermediate track to avoid GC probem
  public FpmEntry inspect4screen(double multiple)
  {
    if (!Tools.isSet(getTrack())) return this;

//    double           ppm0 = 0d, ppm=0d, mz_y1=0d;
//    List<PeakMatch> track = new ArrayList<>(F.getTrack().size());
//    double[]         ppms = new double[F.getTrack().size()], dppm = new double[F.getTrack().size()];
    boolean        has1st = false, y1=false;

    // get a simple average as the starter
//    for (PeakMatch p : F.getTrack()) { ppm0+=p.getMz()*p.getSNR(); ppm+=p.getSNR(); }
//    ppm0/=ppm;

    for (int i=0; i<getTrack().size(); i++)
    {
//      ppm=F.at(i).getMz();
//      PeakMatch pk = new PeakMatch(ppm, Math.abs(F.at(i).getSNR()), F.at(i).getCharge(), ppm-ppm0,F.at(i).getFrequency());
//      pk.setCounts(F.at(i).getCounts());
//      pk.setOriginalMz(F.at(i).getIntensity()).setIndex(i);
//      track.add(pk);
//      ppms[i]=ppm; dppm[i]=(ppm-ppm0); ppm0=ppm;

      // check the presence of y1
      if (at(i).getCharge()==1 && !has1st) has1st=true;
      // the observed m/z was saved in the 'intensity' field!!!
//      if (F.at(i).getCharge()==1 && y1s!=null && y1s.contains(F.at(i).getIntensity()))
//      { y1=true; mz_y1=F.at(i).getIntensity(); }

      // check for the presence of Pro at the N-t of the fragment
      if (at(i).isIonType(IonType.p)) increProlines();
    }
    // capture the individual ppm and remove the worst one if outside 2xsigma
//    double mid = PeakMatch.centroid(track), d0=Stats.mean(dppm, dppm.length),
//        sigma = Stats.stdev(ppms, ppms.length), dsigma=Stats.stdev(dppm, dppm.length);
//
//    int bad=-1, dbad=-1; double worst=multiple*sigma, dworst=multiple*dsigma;
//    for (int i=0; i<track.size(); i++)
//    {
//      if (Math.abs(track.get(i).getMz()-mid)> worst) {  bad=i;  worst=Math.abs(track.get(i).getMz() -mid); }
//      if (Math.abs(track.get(i).getSNR()-d0)>dworst) { dbad=i; dworst=Math.abs(track.get(i).getSNR()-d0); }
//    }

    // return now if there isn't enough matches
//    if (track.size()==0)
//      return F.has1st(has1st).hasExpectedY1(Tools.isSet(y1s) && y1).setIntensity(0).setMotifs(0).setGapScore(0);

    // re-calculate the mid point
//    mid = PeakMatch.centroid(track);

    int best=0, start=0, delta=0, contig_start=0, contig_last=0;
    double scr=0, percentile=0, score=0, sumAI=0d;
    setGapScore0(0);
    for (int i=0; i<getTrack().size(); i++)
    {
      PeakMatch pk = at(i);
//      // correct the drift only if 3 or more fragments are present
//      if (track.size()>2) pk.setMzAndCharge(pk.getMz()-mid, pk.getCharge());
      delta      =  pk.getCharge()-start;
      percentile = (pk.getIntensity()*0.01); // set a minimum
      // accumualte the gap score
      if (delta>0) {
//        scr   =calcGapScore(pk, delta, 1d);
        pk.calcGapScore(delta, 1d);
//        score+=pk.getScore()*percentile;
        score+=pk.getScore() - Math.log10(percentile);
      }
//      else scr=0;

      start=pk.getCharge(); sumAI+=percentile;

      // check the config
      if (i>0 && start-contig_last>1)
      {
        if (contig_last-contig_start>best) best=contig_last-contig_start;
        contig_start=pk.getCharge();
      }
//
//      int first=pk.getCharge(), last=first;
//      for (int j=i+1; j<getTrack().size(); j++)
//        if (last-at(j).getCharge()!=1) break; else last=at(j).getCharge();
//
//      if (first-last>best) best=first-last;
    }

    // http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/tutorials/xegbohtmlnode16.html
    has1st(has1st).setIntensity(sumAI).setMotifs(best).setGapScore(-10d*score);

//    if (track!=null)
//    {
//      for (PeakMatch p : track) p.dispose();
//      track.clear(); track=null;
//    }
//    ppms=dppm=null;

    return this;
  }

  public FpmEntry inspect(Set<Double> y1s, Double local_base)
  {
    if (!Tools.isSet(getTrack())) return this;

    boolean has1st = false, y1=false;
    PeakMatch pk1 = null;
    List<PeakMatch> track = new ArrayList<>(getTrack().size());
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
      PeakMatch pk = track.get(i);

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
    mFragment=(FragmentEntry )Tools.dispose(mFragment);
    mTrack=null;
  }

  @Override
  public void write(DataOutput ds) throws IOException
  {
    IOs.write(ds, mHas1st); IOs.write(ds, mExpectedY1);
    IOs.write(ds, mMotifs); IOs.write(ds, m1stPass); IOs.write(ds, mWeaks); IOs.write(ds, mC13);
    IOs.write(ds, mGapScore); IOs.write(ds, mIntensities); IOs.write(ds, mGapScore0);
//    IOs.write(ds, mMatchScore);
    IOs.write(ds, 0D);
    IOs.write(ds, mFragment);
    IOs.write(ds, mTrack);
  }

  @Override
  public void read(DataInput ds) throws IOException
  {
    mHas1st     =IOs.read(ds,mHas1st);
    mExpectedY1 =IOs.read(ds,mExpectedY1);
    mMotifs     =IOs.read(ds,mMotifs);
    m1stPass    =IOs.read(ds,m1stPass);
    mWeaks      =IOs.read(ds,mWeaks);
    mC13        =IOs.read(ds, mC13);
    mGapScore   =IOs.read(ds, mGapScore);
    mIntensities=IOs.read(ds, mIntensities);
    mGapScore0  =IOs.read(ds, mGapScore0);
//    mMatchScore =IOs.read(ds, mMatchScore);
    IOs.read(ds, 0D); // only a placeholder for compatibility
    mFragment   =IOs.read(ds, mFragment);

    mTrack=IOs.readImmutableList(ds, PeakMatch.class);
  }
}
