package org.ms2ms.data.ms;

import com.google.common.collect.*;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleDirectedWeightedGraph;
import org.ms2ms.algo.Similarity;
import org.ms2ms.data.Binary;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;

/**
 * Created by yuw on 2/24/16.
 */
public class ProteinID implements Comparable<ProteinID>, Binary
{
  public static class DistinctPeptideDesendComparator implements Comparator<ProteinID> {
    public int compare(ProteinID o1, ProteinID o2) {
      return o1 != null && o2 != null ? Double.compare(o2.getSeqMatch().keySet().size(), o1.getSeqMatch().keySet().size()) : 0;
    }
  }

  private Long                           mID;
  private Double                         mBestQVal, mProteoSimilary;
  private String                         mSequence, mAccession, mGene, mName;

  private SimpleDirectedWeightedGraph<SRM, DefaultWeightedEdge> mNetwork;

  private ProteinID                      mParent  =null;
  private Collection<ProteinID>          mChildren=null;
  private Multimap<String, SRMGroup>     mSRMGroups=null;
  private SRMGroup                       mCompositeSRMGroup=null;

  private Multimap<String, PeptideMatch> mSeqMatch=null;
  private Table<String, Integer, PeptideFeature> mSeqChargeFeature = null;

  public ProteinID()                    { super(); }
  public ProteinID(Long s)              { super(); mID=s; }
  public ProteinID(String s, String id) { super(); mSequence=s; mAccession=id; }

  public String getName() { return mName; }
  public String getGene() { return mGene; }
  public Double getBestQVal() { return mBestQVal; }
  public String getSequence() { return mSequence; }
  public String getAccession() { return mAccession; }
  public Double getProteoSimilarity()   { return mProteoSimilary; }

  public Multimap<String, PeptideMatch> getSeqMatch()          { return mSeqMatch; }
  public Multimap<String, SRMGroup>     getSRMGroups()         { return mSRMGroups; }
  public Collection<SRMGroup>           getSRMGroup(String s)  { return mSRMGroups!=null?mSRMGroups.get(s):null; }
  public SRMGroup                       getCompositeSRMGroup() { return mCompositeSRMGroup; }

  public ProteinID setAccession(String s) { mAccession=s; return this; }
  public ProteinID setName(String s) { mName=s; return this; }
  public ProteinID setGene(String s) { mGene=s; return this; }
  public ProteinID setProteoSimilarity(Double s) { mProteoSimilary=s; return this; }
  public ProteinID setCompositeSRMGroup(SRMGroup s) { mCompositeSRMGroup=s; return this; }

  public ProteinID addSRMGroup(SRMGroup g, String s)
  {
    if (mSRMGroups==null) mSRMGroups = HashMultimap.create();
    mSRMGroups.put(s, g);
    return this;
  }
  public int getQualifiedTransitionCounts(float min_apex, float peak_pct)
  {
    int counts=0;
    if (Tools.isSet(getSRMGroups()))
      for (SRMGroup grp : getSRMGroups().values())
        if (Tools.isSet(grp.getSRMs()))
          for (SRM srm : grp.getSRMs().values())
            if (srm.getFragmentMz()>0 && srm.getFeature()!=null && srm.getApex()>=min_apex && srm.getPeakPct()>=peak_pct) counts++;

    return counts;
  }
  public int getExpectedTransitionCounts()
  {
    int counts=0;
    if (Tools.isSet(getSRMGroups()))
      for (SRMGroup grp : getSRMGroups().values())
        if (Tools.isSet(grp.getSRMs()))
          for (SRM srm : grp.getSRMs().values())
            if (srm.getFragmentMz()>0) counts++;

    return counts;
  }

  public ProteinID updateBestQVal(Double s) { if (mBestQVal==null || (s!=null && mBestQVal>s)) mBestQVal=s; return this; }

  public PeptideFeature put(ProteinID pid, PeptideMatch match, String modseq, Integer charge, Double rt, Double mz)
  {
    if (mSeqChargeFeature==null) mSeqChargeFeature = HashBasedTable.create();
    if (mSeqChargeFeature.get(modseq, charge)==null)
    {
      PeptideFeature F = new PeptideFeature(match).addMatch(match).setTitle(modseq).setRT(rt).setMz(mz).setCharge(charge).setProteinID(pid);
      mSeqChargeFeature.put(modseq, charge, F);
    }
    else
    {
      mSeqChargeFeature.get(modseq, charge).addMatch(match);
    }
    return mSeqChargeFeature.get(modseq, charge);
  }
  public PeptideMatch put(String seq, PeptideMatch match)
  {
    if (mSeqMatch==null) mSeqMatch =HashMultimap.create();
    mSeqMatch.put(seq, match);

    return match;
  }
  public PeptideMatch put(PeptideMatch match)
  {
    if (mSeqMatch==null) mSeqMatch =HashMultimap.create();
    mSeqMatch.put(match.toBarePeptide().toString(), match);

    return match;
  }
  public ProteinID putAll(String seq, Collection<PeptideMatch> matches)
  {
    if (seq!=null && Tools.isSet(matches))
      for (PeptideMatch match :matches) put(seq, match);

    return this;
  }

  @Override
  public int hashCode()
  {
    return super.hashCode() + Tools.hashCode(mID,mBestQVal,mSequence,mAccession,mGene) + Tools.hashCode(mSeqMatch);
  }
  @Override
  public String toString()
  {
    return Strs.toString(mID)+"|"+Strs.toString(mAccession)+"|"+Strs.toString(mGene)+" "+Strs.toString(mName)+", "+
        (mSeqMatch!=null?mSeqMatch.size():(mSRMGroups!=null?mSRMGroups.size():0)+
        (Tools.isSet(getSRMGroups())?(", SRM="+getQualifiedTransitionCounts(1000, 75)+"/"+getExpectedTransitionCounts()):"")+
        (getProteoSimilarity()!=null?(", sim="+Tools.d2s(getProteoSimilarity(), 2)):""));
  }

  @Override
  public int compareTo(ProteinID o)
  {
    int r = 0;
    if (o!=null && o.getSequence()!=null && getSequence()!=null)
    {
      r = getSequence().compareTo(o.getSequence());
      if (r==0 && getSeqMatch()!=null && o.getSeqMatch()!=null )
        r = Integer.compare(getSeqMatch().size(), o.getSeqMatch().size());
    }

    return r;
  }
  public static TreeMultimap<Float, ProteinID> toAbundanceRanks(Collection<ProteinID> proteins)
  {
    if (proteins==null) return null;

    // build an index by the key sequences
    TreeMultimap<Float, ProteinID> peptide_protein = TreeMultimap.create(Ordering.natural().reverse(), Ordering.natural().reverse());
    for (ProteinID protein : proteins)
      if (Tools.isSet(protein.getSeqMatch()) && protein.getSequence()!=null)
      {
//        Float abund = (float )protein.getSeqMatch().values().size() / (float )protein.getSequence().length();
        Float abund = (float )protein.getSeqMatch().keySet().size();
        peptide_protein.put(abund, protein);
      }

    return peptide_protein;
  }

  public static Multimap<String, ProteinID> toPeptideProteinIdx(Collection<ProteinID> proteins)
  {
    if (proteins==null) return null;

    // build an index by the key sequences
    Multimap<String, ProteinID> peptide_protein = HashMultimap.create();
    for (ProteinID protein : proteins)
      if (Tools.isSet(protein.getSeqMatch()))
        for (String seq : protein.getSeqMatch().keySet())
          peptide_protein.put(seq, protein);

    return peptide_protein;
  }
  public ProteinID networking() {
    if (getCompositeSRMGroup()==null || !Tools.isSet(getCompositeSRMGroup().getSRMs())) return this;

    // create a new network
    mNetwork = new SimpleDirectedWeightedGraph<>(DefaultWeightedEdge.class);
    List<Float> traces = new ArrayList<>(getCompositeSRMGroup().getSRMs().keySet());

    for (int i = 0; i < traces.size(); i++) {
      SRM x = getCompositeSRMGroup().getSRM(traces.get(i)); mNetwork.addVertex(x);
      for (int j = 0; j < traces.size(); j++) {
        SRM y = getCompositeSRMGroup().getSRM(traces.get(j)); mNetwork.addVertex(y);
        if (i != j) {
          // calc the similarity
          double dp = Similarity.dp_points(x.getXIC(), y.getXIC());
          mNetwork.setEdgeWeight(mNetwork.addEdge(x, y), dp);
        }
      }
    }
    return this;
  }
  @Override
  public void write(DataOutput ds) throws IOException
  {
//    private Long                           mID;
//    private Double                         mBestQVal;
//    private String                         mSequence, mAccession, mGene, mName;
//
//    private ProteinID                      mParent  =null;
//    private Collection<ProteinID>          mChildren=null;
//
//    private Multimap<String, PeptideMatch> mSeqMatch=null;
//    private Table<String, Integer, PeptideFeature> mSeqChargeFeature = null;
  }

  @Override
  public void read(DataInput ds) throws IOException
  {

  }
}
