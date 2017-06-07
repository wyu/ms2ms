package org.ms2ms.data.ms;

import com.google.common.collect.*;
import com.sun.corba.se.impl.encoding.OSFCodeSetRegistry;
import com.thinkaurelius.titan.diskstorage.EntryMetaData;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.ms2ms.algo.PSMs;
import org.ms2ms.utils.Tools;

import java.util.Collection;
import java.util.Map;

/**
 * Created by yuw on 2/24/16.
 */
public class ProteinID implements Comparable<ProteinID>
{
  private Long                           mID;
  private Double                         mBestQVal;
  private String                         mSequence, mAccession, mGene, mName;

  private ProteinID                      mParent  =null;
  private Collection<ProteinID>          mChildren=null;

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
  public Multimap<String, PeptideMatch> getSeqMatch() { return mSeqMatch; }

  public ProteinID setAccession(String s) { mAccession=s; return this; }
  public ProteinID setName(String s) { mName=s; return this; }
  public ProteinID setGene(String s) { mGene=s; return this; }

  public ProteinID updateBestQVal(Double s) { if (mBestQVal==null || (s!=null && mBestQVal>s)) mBestQVal=s; return this; }

  public PeptideFeature put(PeptideMatch match, String modseq, Integer charge)
  {
    if (mSeqChargeFeature==null) mSeqChargeFeature = HashBasedTable.create();
    if (mSeqChargeFeature.get(modseq, charge)==null)
    {
      PeptideFeature F = new PeptideFeature(match).addMatch(match).setTitle(modseq);
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

    // build an index by the peptide sequences
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

    // build an index by the peptide sequences
    Multimap<String, ProteinID> peptide_protein = HashMultimap.create();
    for (ProteinID protein : proteins)
      if (Tools.isSet(protein.getSeqMatch()))
        for (String seq : protein.getSeqMatch().keySet())
          peptide_protein.put(seq, protein);

    return peptide_protein;
  }
}
