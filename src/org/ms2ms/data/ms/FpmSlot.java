package org.ms2ms.data.ms;

import com.google.common.collect.ImmutableList;
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
public class FpmSlot implements Comparable<FpmSlot>, Disposable, Binary
{
  private char                 mMotifs=0;
  private float              mGapScore=0;
//  private Double              mMatchScore =null;
  private FragmentSlot       mFragment   =null;
  private ImmutableList<PeakEntry> mTrack =null;

  public static class PositionComparator implements Comparator<FpmSlot>
  {
    public PositionComparator() { }
    public int compare(FpmSlot o1, FpmSlot o2)
    {
      return Integer.compare(Math.abs(o1.getFragment().getPeptideKey()), Math.abs(o2.getFragment().getPeptideKey()));
    }
  }
  public static class GapScoreDecendComparator implements Comparator<FpmSlot>
  {
    public GapScoreDecendComparator() { }
    public int compare(FpmSlot o1, FpmSlot o2)
    {
      return Double.compare(o2.getGapScore(), o1.getGapScore());
    }
  }

  public FpmSlot()
  {
    super(); mFragment = new FragmentSlot();
  }
  public FpmSlot(FpmSlot f)
  {
    super();
    if (f!=null)
    {
      mMotifs=f.mMotifs;
      mGapScore=f.mGapScore;
      mFragment = new FragmentSlot(f.mFragment);
      mTrack = f.mTrack!=null?f.mTrack.asList():null;
    }
  }
  public FpmSlot(FragmentSlot f, List<PeakEntry> t)
  {
    super();
    mFragment=f; mTrack=ImmutableList.copyOf(t);
  }
  public FpmSlot(FragmentSlot f, PeakEntry[] matches, int size)
  {
    super();
    setTrack(matches, 0, size);
    mFragment=f;
  }

  public boolean             isDecoy()      { return mFragment!=null && mFragment.getPeptideKey()<0; }
  public double              getGapScore()  { return mGapScore; }
  public int                 getMotifs()    { return mMotifs; }
  public FragmentSlot        getFragment()  { return mFragment; }
  public ImmutableList<PeakEntry> getTrack()     { return mTrack; }
  public PeakEntry       at(int s)      { return mTrack.get(s); }

  public FpmSlot setMotifs(int s) { mMotifs = (char )s;  return this; }

  public FpmSlot setGapScore(double s) { mGapScore  = (float)s; return this; }
  public FpmSlot setTrack(ImmutableList<PeakEntry> s) { mTrack=s; return this; }

  public FpmSlot setTrack(PeakEntry[] track, int beg, int end)
  {
    List<PeakEntry> tr = new ArrayList<>(end-beg);
    for (int i=beg; i<end; i++) tr.add(new PeakEntry(track[i]));

    mTrack = ImmutableList.copyOf(tr);
    return this;
  }
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
  public FpmSlot shallow_copy()
  {
    FpmSlot clone = new FpmSlot();
    // no scores, just the skeleton
    clone.mFragment = mFragment;
    clone.mTrack=mTrack;

    return clone;
  }

//  // update the ion matches and return the local base intensity
//  public double match(OffsetPpmTolerance tol)
//  {
//    double local_base=0;
//    if (Tools.isSet(getTrack()))
//    {
//      List<PeakMatch> pts = new ArrayList<>(getTrack().size());
//      for (PeakMatch  pt : getTrack())
//      {
//        Double calc = pt.getScore(), obs = pt.getIntensity();
//        if (calc!=null && obs!=null && tol.withinTolerance(calc, obs))
//        {
//          PeakMatch p = pt.clone();
//          p.setMzAndCharge(Stats.ppm(obs, calc)+tol.getOffset(obs), pt.getCharge());
//          pts.add(p);
//          if (p.getSNR()>local_base) local_base=p.getSNR();
//        }
//      }
//      mTrack = ImmutableList.copyOf(pts);
//    }
//
//    return local_base;
//  }
//  // following the re-match, update the scores of the series
//  public FpmSlot score(Double local_base)
//  {
//    if (Tools.isSet(getTrack()))
//    {
//      // we used to score only the y1.
//      int best=0, start=1, delta=0;
//      double percentile=0, score=(has1st()?calcGapScore(getTrack().get(getTrack().size()-1), 1, 1d)*0.1d:0),
//          sumAI=0d, base=local_base!=null?1d/local_base:0.01d;
//      for (int i=getTrack().size()-1; i>=0; i--)
//      {
//        PeakMatch pk = getTrack().get(i);
//
//        delta = pk.getCharge()-start;
//        percentile = (pk.getSNR()*base); // set a minimum
//        // accumualte the gap score
//        if (delta>0) score+=(calcGapScore(pk, delta, 1d))*percentile;
//
//        start=pk.getCharge(); sumAI+=percentile;
//
//        int first=getTrack().get(i).getCharge(), last=first;
//        for (int j=i+1; j<getTrack().size(); j++)
//          if (last-getTrack().get(j).getCharge()!=1) break; else last=getTrack().get(j).getCharge();
//
//        if (first-last>best) best=first-last;
//      }
//      setMotifs(best).setGapScore(-10d*score);
//    }
//
//    return this;
//  }
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
  public int compareTo(FpmSlot o)
  {
    int c = mFragment.compareTo(o.getFragment());
    if (c==0 && mTrack!=null && o.getTrack()!=null) c = Integer.compare(mTrack.size(), o.getTrack().size());

    return c;
  }

  @Override
  public String toString()
  {
    return "#"+(Tools.isSet(mTrack)?getMotifs()+"/"+mTrack.size():0)+"@"+(mFragment!=null?mFragment.toString():"")+"$"+Tools.d2s(getGapScore(), 1);
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
//  public Collection<Double> fillMz(Collection<Double> data)
//  {
//    if (Tools.isSet(getTrack()))
//      for (PeakMatch p : getTrack()) data.add(p.getIntensity());
//
//    return data;
//  }
//  public FpmSlot inspect(Set<Double> y1s, Double local_base)
//  {
//    if (!Tools.isSet(getTrack())) return this;
//
//    boolean has1st = false, y1=false;
//    PeakMatch pk1 = null;
//    List<PeakMatch> track = new ArrayList<>(getTrack().size());
//    for (int i=0; i<getTrack().size(); i++)
//    {
//      track.add(at(i));
//      // check the presence of y1
//      if (at(i).getCharge()==1 && !has1st) { has1st=true; pk1=at(i); }
//      // the observed m/z was saved in the 'intensity' field!!!
//      if (at(i).getCharge()==1 && y1s!=null && y1s.contains(at(i).getIntensity())) y1=true;
//    }
//
//    // return now if there isn't enough matches
//    if (track.size()==0)
//      return has1st(has1st).hasExpectedY1(Tools.isSet(y1s) && y1).setIntensity(0).setMotifs(0).setGapScore(0);
//
//    // we used to score only the y1.
//    int best=0, start=1, delta=0;
//    double percentile=0, score=(has1st?calcGapScore(track.get(track.size() - 1), 1, 1d)*0.1d:0),
//        sumAI=0d, base=local_base!=null?1d/local_base:0.01d;
//    for (int i=track.size()-1; i>=0; i--)
//    {
//      PeakMatch pk = track.get(i);
//
//      delta = pk.getCharge()-start;
//      percentile = (pk.getSNR()*base); // set a minimum
//      // accumualte the gap score
//      if (delta>0) score+=(calcGapScore(pk, delta, 1d))*percentile;
//
//      start=pk.getCharge(); sumAI+=percentile;
//
//      int first=track.get(i).getCharge(), last=first;
//      for (int j=i+1; j<track.size(); j++)
//        if (last-track.get(j).getCharge()!=1) break; else last=track.get(j).getCharge();
//
//      if (first-last>best) best=first-last;
//    }
//
//    // http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/tutorials/xegbohtmlnode16.html
//    return has1st(has1st).hasExpectedY1(Tools.isSet(y1s) && y1).setIntensity(sumAI).setMotifs(best).setGapScore(-10d*score);
//  }

  @Override
  public FpmSlot clone()
  {
    return new FpmSlot(this);
  }
  @Override
  public int hashCode()
  {
    int hash = mMotifs+Tools.hashCodes(mGapScore);
    if (mFragment==null) hash+=mFragment.hashCode();

    return hash+hashcodeByTrack();
  }

  @Override
  public void dispose()
  {
    mFragment=(FragmentSlot )Tools.dispose(mFragment);
    mTrack=null;
  }

  @Override
  public void write(DataOutput ds) throws IOException
  {
    IOs.write(ds, mMotifs);
    IOs.write(ds, mGapScore);
    IOs.write(ds, mFragment);
    IOs.write(ds, mTrack);
  }

  @Override
  public void read(DataInput ds) throws IOException
  {
    mMotifs     =IOs.read(ds,mMotifs);
    mGapScore   =IOs.read(ds, mGapScore);
    mFragment   =IOs.read(ds, mFragment);

    mTrack=IOs.readImmutableList(ds, PeakEntry.class);
  }
}
