package org.ms2ms.data.ms;

import com.google.common.collect.*;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.ms2ms.data.Binary;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.Collection;
import java.util.Comparator;

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

  private ProteinID                      mParent  =null;
  private Collection<ProteinID>          mChildren=null;
  private Multimap<String, SRMGroup>     mSRMGroups=null;

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

  public Multimap<String, PeptideMatch> getSeqMatch()         { return mSeqMatch; }
  public Multimap<String, SRMGroup>     getSRMGroups()        { return mSRMGroups; }
  public Collection<SRMGroup>           getSRMGroup(String s) { return mSRMGroups!=null?mSRMGroups.get(s):null; }

  public ProteinID setAccession(String s) { mAccession=s; return this; }
  public ProteinID setName(String s) { mName=s; return this; }
  public ProteinID setGene(String s) { mGene=s; return this; }
  public ProteinID setProteoSimilarity(Double s) { mProteoSimilary=s; return this; }

  public ProteinID addSRMGroup(SRMGroup g, String s)
  {
    if (mSRMGroups==null) mSRMGroups = HashMultimap.create();
    mSRMGroups.put(s, g);
    return this;
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
    return Strs.toString(mID)+"|"+Strs.toString(mAccession)+"|"+Strs.toString(mGene)+" "+Strs.toString(mName)+", "+(mSeqMatch!=null?mSeqMatch.size():0);
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
