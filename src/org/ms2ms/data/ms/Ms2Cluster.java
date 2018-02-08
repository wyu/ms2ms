package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.algo.Peaks;
import org.ms2ms.algo.Similarity;
import org.ms2ms.algo.Spectra;
import org.ms2ms.io.BufferedRandomAccessFile;
import org.ms2ms.io.MsIO;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Tools;

import java.io.IOException;
import java.util.*;

public class Ms2Cluster implements Comparable<Ms2Cluster>
{
  private float mMz, mRT;
  private int mByMz=0, mByMzRT=0, mByMzRtFrag=0, mCharge=0;
  private String mName;

  private MsnSpectrum mMaster; // a composite spectrum to represent the cluster

  private Ms2Pointer mHead;
  private Collection<Ms2Pointer> mMembers = new TreeSet<>(), mCandidates = new TreeSet<>(); // actual or possible members of the cluster


  public Ms2Cluster() { super(); }
  public Ms2Cluster(String s) { super(); mName=s; }
  public Ms2Cluster(Ms2Pointer s) { super(); mHead=s; }

  public Collection<Ms2Pointer> getCandidates() { return mCandidates; }
  public Collection<Ms2Pointer> getMembers()    { return mMembers; }
  public MsnSpectrum            getMaster()     { return mMaster; }
  public Ms2Pointer             getHead()       { return mHead; }

  public int    getCandidateSize() { return mCandidates!=null?mCandidates.size():0; }; // the head is already a part of the candidates and members
  public int    size()           { return mMembers!=null?mMembers.size():0; }; // the head is already a part of the candidates and members
  public int    getNbyMz()       { return mByMz; }
  public int    getNbyMzRT()     { return mByMzRT; }
  public int    getCharge()      { return mCharge; }
  public float  getMz()          { return mMz; }
  public float  getRT()          { return mRT; }
  public int    getNbyMzRtFrag() { return mByMzRtFrag; }
  public String getName()        { return mName; }

  public Ms2Cluster setHead(Ms2Pointer s) { mHead      =s; return this; }
  public Ms2Cluster setNbyMz(      int s) { mByMz      =s; return this; }
  public Ms2Cluster setNbyMzRt(    int s) { mByMzRT    =s; return this; }
  public Ms2Cluster setNbyMzRtFrag(int s) { mByMzRtFrag=s; return this; }
  public Ms2Cluster setName(    String s) { mName      =s; return this; }
  public Ms2Cluster setMz(float s) { mMz=s; return this; }
  public Ms2Cluster setRT(float s) { mRT=s; return this; }

  public Ms2Cluster addCandidate(Ms2Pointer s) { if (s!=null) mCandidates.add(s); return this; }
  public Ms2Cluster addMember(   Ms2Pointer s) { if (s!=null) mMembers.add(s); return this; }

  public int getCandidateRemain() { return (mCandidates!=null?mCandidates.size():0)-(mMembers!=null?mMembers.size():0); }
  public boolean contains(Ms2Pointer p)
  {
    if (Tools.equals(p, mHead) || Tools.contains(mCandidates, p)) return true;
    return false;
  }
  public Ms2Cluster setMaster(BufferedRandomAccessFile bin) throws IOException
  {
    if (size()<=1) mMaster = MsIO.readSpectrumIdentifier(bin, new MsnSpectrum(), getHead().pointer);
    return this;
  }
  public Ms2Cluster resetMembers()
  {
    if (mMembers!=null) mMembers.clear(); else mMembers = new TreeSet<>();
    return this;
  }
  public Ms2Cluster resetCandidates()
  {
    if (mCandidates!=null) mCandidates.clear(); else mCandidates = new TreeSet<>();
    return this;
  }
  /// start the method section  ///

  // from the HEAD to the candidates at high cutoff
  public Ms2Cluster cluster(Map<Ms2Pointer, MsnSpectrum> spectra, Float lowmass, OffsetPpmTolerance tol, double min_dp)
  {
    // quit if the master or input are not present!
    if (mHead==null || !Tools.isSet(spectra)) return null;

    // setting the seed spectrum
    MsnSpectrum HEAD = (mMaster!=null?mMaster:spectra.get(mHead));

    // the collections
    if (mMembers!=null) mMembers.clear(); else mMembers = new ArrayList<>();
    Collection<MsnSpectrum> members = new ArrayList<>();
    // setup the master first
    List<Peak> head = Spectra.toListOfPeaks(HEAD, lowmass);
    for (Ms2Pointer member : mCandidates)
    {
      MsnSpectrum scan = spectra.get(member);
      if (scan!=null)
      {
        // calc the forward and backward DPs and choose the smallest
        member.dp=(float )Similarity.bidirectional_dp(head, Spectra.toListOfPeaks(scan,lowmass), tol, true, true, true);
        if (member.dp>=min_dp)
        {
          mMembers.add(member);
          members.add(spectra.get(member));
        }
      }
    }
    if (members.size()>1)
    {
      mMaster = Spectra.accumulate(HEAD, tol, 0.5f, members);
      mMz     = (float )mMaster.getPrecursor().getMz();
      mHead.cluster=this;
    }

    // remove the local objects
    head = (List )Tools.dispose(head);
    members = Tools.dispose(members);

    return this;
  }
  public Ms2Cluster calcMz()
  {
    if (Tools.isSet(mMembers))
    {
      Collection<Double> ms = new ArrayList<>(); double z=0;
      for (Ms2Pointer p : mMembers) { ms.add((double )p.mz); z+=p.z; }
      mMz = (float )Stats.mean(ms); mCharge=(int )Math.round(z);
    }
    return this;
  }
  // trim away the members in the cluster by the matching probability to the 'center'.
  public Ms2Cluster trimByMatchProb(Map<Ms2Pointer, MsnSpectrum> spectra, int regions, double cut)
  {
    // quit if the master or input are not present!
    if (mMaster==null || !Tools.isSet(spectra)) return null;

    // grab the index from our master first
    List<Peak> master = Similarity.index(Spectra.toListOfPeaks(mMaster), regions, 1, 5, 0);
    long bins = (long )(Math.log(2)/Math.log(1d+5E-6)), npks = master.size();

    Iterator<Ms2Pointer> itr = mCandidates.iterator();
    while (itr.hasNext())
    {
      MsnSpectrum scan = spectra.get(itr.next());
      if (scan!=null)
      {
        List<Peak> member = Similarity.index(Spectra.toListOfPeaks(scan), regions, 1, 5, 0);
        int overlap = Peaks.overlap_counts(master, member, 0.01, true);

        double prob = Stats.hypergeom(overlap, npks, member.size(), bins);

        // check the exit condition
        if (-1d*prob<cut) itr.remove();
        // clean up the objects
        member = (List )Tools.dispose(member);
      }
    }

    return this;
  }
  // trim the matching candidates by matching probability
  public Ms2Cluster trimByDotP(Map<Ms2Pointer, MsnSpectrum> spectra, Tolerance tol, double min_dp)
  {
    // quit if the master or input are not present!
    if (mMaster==null || !Tools.isSet(spectra)) return null;

    // setup the master first
    List<Peak>        master = Spectra.toListOfPeaks(mMaster);
    Iterator<Ms2Pointer> itr = mCandidates.iterator();
    while (itr.hasNext())
    {
      MsnSpectrum scan = spectra.get(itr.next());
      if (scan!=null && Similarity.dp(master, Spectra.toListOfPeaks(scan), tol, true, true)<min_dp) itr.remove();
    }

    return this;
  }
  @Override
  public int hashCode()
  {
    int hc = mHead!=null?mHead.hashCode():0;
    if (mName!=null) hc+=mName.hashCode();

    hc += mByMz+mByMzRT+mByMzRtFrag;

    if (Tools.isSet(mCandidates))
      for (Ms2Pointer p : mCandidates) hc+=p.hashCode();

    if (Tools.isSet(mMembers))
      for (Ms2Pointer p : mMembers) hc+=p.hashCode();

    return hc;
  }

  @Override
  public int compareTo(Ms2Cluster o)
  {
    return Integer.compare(hashCode(), o.hashCode());
  }
}
