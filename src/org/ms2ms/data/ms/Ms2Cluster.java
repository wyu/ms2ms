package org.ms2ms.data.ms;

import com.google.common.collect.*;
import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.Disposable;
import org.ms2ms.algo.Peaks;
import org.ms2ms.algo.PurgingPeakProcessor;
import org.ms2ms.algo.Similarity;
import org.ms2ms.algo.Spectra;
import org.ms2ms.data.Binary;
import org.ms2ms.r.Dataframe;
import org.ms2ms.utils.BufferedRandomAccessFile;
import org.ms2ms.io.MsIO;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.*;

public class Ms2Cluster implements Comparable<Ms2Cluster>, Binary, Disposable, Ion
{
  @Override
  public void dispose()
  {
    mMaster     = null;
    mHead       = null;
    mMembers    = Tools.dispose(mMembers);
    mCandidates = Tools.dispose(mCandidates);
  }

  public enum NodeType { REF, MSMS, CELL, NULL, NONE };

  private NodeType mType = NodeType.MSMS;
  private float mMz, mRT, mImpurity=-1;
  private int mByMz=0, mByMzRT=0, mByMzRtFrag=0, mCharge=0, mNamed=0;
  private String mName="", mID="", mMajority=null;

  private SortedSetMultimap<Float, String> mPctName;

  private MsnSpectrum mMaster; // a composite spectrum to represent the cluster

  private Ms2Pointer mHead, mMasterHead;
  private Collection<Ms2Pointer> mMembers = new TreeSet<>(), mCandidates = new TreeSet<>(); // actual or possible members of the cluster

  // no need to save to the archive
  private SortedMap<Float,Float> mMasterIonMap;
  private Set<Float> mIndexIons;

  public Ms2Cluster() { super(); }
  public Ms2Cluster(String s) { super(); mName=s; }
  public Ms2Cluster(String n, String id, NodeType t) { super(); mName=n; mID=id; mType=t; }
  public Ms2Cluster(Ms2Pointer s) { super(); mHead=s; }

  public Collection<Ms2Pointer> getCandidates() { return mCandidates; }
  public Collection<Ms2Pointer> getMembers()    { return mMembers; }
  public MsnSpectrum            getMaster()     { return mMaster; }
  public Ms2Pointer             getHead()       { return mHead; }
  public Ms2Pointer             getMasterHead() { return mMasterHead; }

  public boolean isType(NodeType... t) { return Tools.isA(mType, t); }
  public boolean contains(String run, Integer scans)
  {
    if (Tools.isSet(getCandidates()))
      for (Ms2Pointer p : getCandidates())
        if (p.run.equals(run) &&
          Tools.isA(p.scan,scans)) return true;

    return false;
  }
  public int    getCandidateSize() { return mCandidates!=null?mCandidates.size():0; }; // the head is already a part of the candidates and members
  public int    size()           { return mMembers!=null?mMembers.size():0; }; // the head is already a part of the candidates and members
  public int    getNbyMz()       { return mByMz; }
  public int    getNbyMzRT()     { return mByMzRT; }
  public int    getCharge()      { return mCharge; }
  public float  getMz()          { return mMz; }
  public float  getMH()          { return Peaks.toMH(getMz(),getCharge()); }
  public float  getRT()          { return mRT; }
  public float  getImpurity()    { return mImpurity; }
  public int    getNbyMzRtFrag() { return mByMzRtFrag; }
  public String getName()        { return mName; }
  public String getID()          { return mID; }
  public String getMajority()    { return mMajority; }
  public NodeType getType()      { return mType; }

  public Ms2Cluster setType(NodeType   s) { mType      =s; return this; }
  public Ms2Cluster setHead(Ms2Pointer s) { mHead      =s; return this; }
  public Ms2Cluster setCharge(     int s) { mCharge    =s; return this; }
  public Ms2Cluster setNbyMz(      int s) { mByMz      =s; return this; }
  public Ms2Cluster setNbyMzRt(    int s) { mByMzRT    =s; return this; }
  public Ms2Cluster setNbyMzRtFrag(int s) { mByMzRtFrag=s; return this; }
  public Ms2Cluster setName(    String s) { mName      =s; return this; }
  public Ms2Cluster setID(      String s) { mID        =s; return this; }
  public Ms2Cluster setMajorityID(String s) { mMajority=s; return this; }
  public Ms2Cluster setMz(       float s) { mMz        =s; return this; }
  public Ms2Cluster setRT(       float s) { mRT        =s; return this; }
  public Ms2Cluster setImpurity( float s) { mImpurity  =s; return this; }
  public Ms2Cluster setMasterHead(Ms2Pointer s) { mMasterHead=s; return this; }

  public Ms2Cluster addCandidate(Ms2Pointer s) { if (s!=null) mCandidates.add(s); return this; }
  public Ms2Cluster addMember(   Ms2Pointer s) { if (s!=null) mMembers.add(s); return this; }
  public Ms2Cluster addMembers(Collection<Ms2Pointer> s) { if (s!=null) mMembers.addAll(s); return this; }

  public Ms2Cluster updateMajorityByName()
  {
    if (Tools.isSet(getMembers()))
    {
      Multimap<String, Ms2Pointer> name_mem = HashMultimap.create();
      for (Ms2Pointer p : getMembers())
        Tools.put(name_mem, Tools.back(Strs.split(p.name,'$')), p);

      if (mPctName==null) mPctName = TreeMultimap.create(Ordering.natural().reverse(), Ordering.natural()); else mPctName.clear();
      for (String name : name_mem.keySet())
        mPctName.put(name_mem.get(name).size()*100f/size(), name);

      mMajority=null; mImpurity = Collections.max(mPctName.keySet());
      for (Float pct : mPctName.keySet())
        for (String name : mPctName.get(pct))
          mMajority = Strs.extend(mMajority, Tools.d2s(pct,1)+"%"+name,";");

      name_mem = Tools.dispose(name_mem);
    }
    return this;
  }
  public int getNamed()
  {
    mNamed=0;
    if (Tools.isSet(getCandidates()))
      for (Ms2Pointer p : getCandidates())
        if (Strs.isSet(p.name)) mNamed++;

    return mNamed;
  }
  public Set<Float> indexFromIonMap()
  {
    mIndexIons = Similarity.index(mMasterIonMap, 7, 0,1,5,0);
    return mIndexIons;
  }
  public SortedMap<Float,Float> getMasterIonMap() { return mMasterIonMap; }
  public Set<Float> getIndexIons() { return mIndexIons; }
  public int getCandidateRemain() { return (mCandidates!=null?mCandidates.size():0)-(mMembers!=null?mMembers.size():0); }
  public boolean contains(Ms2Pointer p)
  {
    if (Tools.equals(p, mHead) || Tools.contains(mCandidates, p)) return true;
    return false;
  }
  public Ms2Cluster setMaster(MsnSpectrum s) { mMaster=s; return this; }
  public Ms2Cluster setMaster(BufferedRandomAccessFile bin) throws IOException
  {
    if (size()<=1) mMaster = MsIO.readSpectrumIdentifier(bin, new MsnSpectrum(), getHead().pointer);
    return this;
  }
  public Ms2Cluster prepMaster(BufferedRandomAccessFile bin) throws IOException
  {
    if (size()>0)
    {
      mMaster = new MsnSpectrum(); // reset the master
      for (Ms2Pointer p : getMembers())
        mMaster.addPeaks(MsIO.readSpectrumIdentifier(bin, new MsnSpectrum(), p.pointer));
    }

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
  public Ms2Cluster cluster(Map<Ms2Pointer, MsnSpectrum> spectra, Float lowmass, OffsetPpmTolerance tol, double min_dp, int miss_index)
  {
    // quit if the master or input are not present!
    if (mHead==null || !Tools.isSet(spectra)) return null;

    mHead.cluster=this;

    // setting the seed spectrum
    if (getMaster()==null && spectra.get(getHead())!=null)
      setMaster(spectra.get(getHead()).copy(new PurgingPeakProcessor()));

    // setup the master first
    List<Peak> head = Spectra.toListOfPeaks(getMaster(), lowmass);
    // get the index of the HEAD spectrum
    List<Peak> index = Similarity.index(head, 7, 0, 1, 5, 0);

    // the collections
    if (mMembers!=null) mMembers.clear(); else mMembers = new ArrayList<>();
    // new members from this round
    Collection<MsnSpectrum> members = new ArrayList<>();
    // the error tolerance
    double delta = tol.calcError(500);
    for (Ms2Pointer member : mCandidates)
    {
      MsnSpectrum scan = spectra.get(member);
      // need to exclude the head itself
      if (scan!=null && member.hcode!=getHead().hcode)
      {
        // calc the forward and backward DPs and choose the smallest
        List<Peak> pks = Spectra.toListOfPeaks(scan,lowmass);
        // make sure a min number of the index peaks are found
        if (miss_index>=0 && index.size()-Peaks.overlap_counts(pks, index, delta, true)>miss_index) continue;

        float dp = (float )Similarity.bidirectional_dp(head, pks, tol, true, true, true);

        if (dp>=min_dp)
        {
          member.dp=dp;
          // now the matching probability
          member.prob = (float )Similarity.similarity_hg(head, pks, delta);

          member.cluster=this;
          mMembers.add(member);
          members.add(spectra.get(member));
        }
      }
    }
    mMaster = Spectra.accumulate(getMaster(), tol, 0.5f, members);
    mMz     = (float )mMaster.getPrecursor().getMz();
    mRT     = (float )getMaster().getRetentionTimes().getFirst().getTime()/60f;
    mCharge = getMaster().getPrecursor().getCharge();
//    mID     = Tools.d2s(getMz(),3)+"|"+getCharge()+"|"+Tools.d2s(getRT(),1);
    mID     = toString();

    // remove the local objects
    head    = (List )Tools.dispose(head);
    index   = (List )Tools.dispose(index);
    members = Tools.dispose(members);

    return this;
  }
  // for exact match without mass tolerance
  public Ms2Cluster cluster(Map<Ms2Pointer, Map<Float,Float>> spectra, double min_dp)
  {
    // quit if the master or input are not present!
    if (mHead==null || !Tools.isSet(spectra)) return null;

    mHead.cluster=this;

    Map<Float,Float> head = spectra.get(getHead());

    if (mMasterIonMap==null) mMasterIonMap = new TreeMap<>(head);

    // the collections
    if (mMembers!=null) mMembers.clear(); else mMembers = new ArrayList<>();
    for (Ms2Pointer member : mCandidates)
    {
      Map<Float,Float> scan = spectra.get(member);
      // need to exclude the head itself
      if (scan!=null && member.hcode!=getHead().hcode)
      {
        float dp = (float )Similarity.bidirectional_dp(head, scan, true,true, min_dp);
        if (dp>=min_dp)
        {
          member.dp=dp;
          // now the matching probability
          member.prob = (float )Similarity.similarity_hg(head, scan, 25000);

          member.cluster=this;
          mMembers.add(member);
          mMasterIonMap = (SortedMap )Tools.accumulate(mMasterIonMap, scan);
        }
      }
    }
    mID = toString();

    return this;
  }

  public Ms2Cluster calcMz()
  {
    if (Tools.isSet(mMembers))
    {
      Collection<Double> ms = new ArrayList<>(); double z=0;
      for (Ms2Pointer p : mMembers) { ms.add((double )p.mz); z+=p.z; }
      mMz = (float )Stats.mean(ms);
      mCharge=(int )Math.round((double )z/(double )getMembers().size());
    }
    return this;
  }
  // trim away the members in the cluster by the matching probability to the 'center'.
  public Ms2Cluster trimByMatchProb(Map<Ms2Pointer, MsnSpectrum> spectra, int regions, double cut)
  {
    // quit if the master or input are not present!
    if (mMaster==null || !Tools.isSet(spectra)) return null;

    // grab the index from our master first
    List<Peak> master = Similarity.index(Spectra.toListOfPeaks(mMaster), regions, 0, 1, 5, 0);
    long bins = (long )(Math.log(2)/Math.log(1d+5E-6)), npks = master.size();

    Iterator<Ms2Pointer> itr = mCandidates.iterator();
    while (itr.hasNext())
    {
      MsnSpectrum scan = spectra.get(itr.next());
      if (scan!=null)
      {
        List<Peak> member = Similarity.index(Spectra.toListOfPeaks(scan), regions, 0, 1, 5, 0);
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
  public void updateAnnotations(Map<String, Ms2Cluster> runscan_pepcls)
  {
    if (size()>15)
      System.out.println();
    if (Tools.isSet(getCandidates()) && Tools.isSet(runscan_pepcls))
      for (Ms2Pointer p : getCandidates())
      {
        String runscan = p.run+"#"+p.scan;
        if (runscan_pepcls.get(runscan)!=null)
          p.name = runscan_pepcls.get(runscan).getName();
      }
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
  public Dataframe details(Dataframe df, String rowid)
  {
    df.put(rowid, "Majority", mMajority);
    df.put(rowid, "ID", mID);
    if (getMaster()!=null)
    {
      if (getMaster().getPrecursor()!=null)
      {
        df.put(rowid, "PrecMH", Tools.d2s(Peaks.toMH(getMaster().getPrecursor()), 4));
        df.put(rowid, "PrecZ",  getMaster().getPrecursor().getCharge());
      }
      if (Tools.isSet(getMaster().getRetentionTimes()) && getMaster().getRetentionTimes().getFirst()!=null)
        df.put(rowid, "RT", Tools.d2s(getMaster().getRetentionTimes().getFirst().getTime()/60d, 2));
    }
    if (getMembers()   !=null) df.put(rowid, "mN", getMembers().size()+"");
    if (getCandidates()!=null) df.put(rowid, "cN", getCandidates().size()+"");
    if (Strs.isSet(getName())) df.put(rowid, "Name", getName());

    if (Tools.isSet(mPctName))
    {
      int n=0;
      for (Float pct : mPctName.keySet())
        for (String name : mPctName.get(pct))
        {
          if (++n>3) break;
          df.put(rowid, "Pct."+n, Tools.d2s(pct, 1));
          df.put(rowid, "Tag."+n, name);
        }
    }
    df.put(rowid, "Impurity", Tools.d2s(mImpurity,2));
    return df;
  }
  @Override
  public String toString()
  {
    String s=mMajority;
    if (getMaster()!=null)
    {
      s+=";";
      if (getMaster().getPrecursor()!=null)
        s = Strs.extend(s, (Tools.d2s(Peaks.toMH(getMaster().getPrecursor()), 4)+"|"+ getMaster().getPrecursor().getCharge()+"|"), "");
      if (Tools.isSet(getMaster().getRetentionTimes()) && getMaster().getRetentionTimes().getFirst()!=null)
        s = Strs.extend(s, Tools.d2s(getMaster().getRetentionTimes().getFirst().getTime()/60d, 2), "");
    }
    s = Strs.extend(s,
        (getMembers()!=null?getMembers().size():0) + "/"+
        (getCandidates()!=null?getCandidates().size():0) +
        (Strs.isSet(getName())?"::"+getName():"") + "|"+hashCode(), "");

    if (!Strs.isSet(mMajority))
      s = Strs.extend(s,
        (Strs.isSet(getName())?"::"+getName():"") + "|"+hashCode(), "");

    return s;
  }
  @Override
  public int compareTo(Ms2Cluster o)
  {
    return Integer.compare(hashCode(), o.hashCode());
  }

  @Override
  public void write(DataOutput ds) throws IOException
  {
    IOs.write(ds, mType.name());
    IOs.write(ds, mMz);
    IOs.write(ds, mRT);
    IOs.write(ds, mImpurity);
    IOs.write(ds, mByMz);
    IOs.write(ds, mByMzRT);
    IOs.write(ds, mByMzRtFrag);
    IOs.write(ds, mCharge);
    IOs.write(ds, mName);
    IOs.write(ds, mID);
    IOs.write(ds, mMajority);
    IOs.write(ds, mHead);

    IOs.write(ds, mMembers);
    IOs.write(ds, mCandidates);

    // write the size of the peaks
    IOs.write(ds, mMaster!=null?mMaster.size():0);
    if (mMaster!=null) MsIO.write(ds, mMaster);
  }

  @Override
  public void read(DataInput ds) throws IOException
  {
    mType       = NodeType.valueOf(IOs.read(ds, ""));
    mMz         = IOs.read(ds, 0f);
    mRT         = IOs.read(ds, 0f);
    mImpurity   = IOs.read(ds, 0f);
    mByMz       = IOs.read(ds, 0);
    mByMzRT     = IOs.read(ds, 0);
    mByMzRtFrag = IOs.read(ds, 0);
    mCharge     = IOs.read(ds, 0);

    mName       = IOs.read(ds, "");
    mID         = IOs.read(ds, "");
    mMajority   = IOs.read(ds, "");
    mHead       = IOs.read(ds, new Ms2Pointer());

    mMembers    = IOs.readList(ds, Ms2Pointer.class);
    mCandidates = IOs.readList(ds, Ms2Pointer.class);

    int npks = IOs.read(ds, 0);
    if (npks>0) mMaster = MsIO.readSpectrumIdentifier(ds, new MsnSpectrum());

    if (!Strs.isSet(mID)) mID = toString();
  }
}
