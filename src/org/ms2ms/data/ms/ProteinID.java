package org.ms2ms.data.ms;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Ordering;
import com.google.common.collect.TreeMultimap;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.ms2ms.utils.Tools;

import java.util.Collection;

/**
 * Created by yuw on 2/24/16.
 */
public class ProteinID implements Comparable<ProteinID>
{
  private Double                         mBestQVal;
  private String                         mSequence, mAccessionId;

  private ProteinID                      mParent  =null;
  private Collection<ProteinID>          mChildren=null;

  private Multimap<String, PeptideMatch> mSeqMatch=null;

  public ProteinID() { super(); }
  public ProteinID(String s, String id) { super(); mSequence=s; mAccessionId=id; }

  public Double getBestQVal() { return mBestQVal; }
  public String getSequence() { return mSequence; }
  public Multimap<String, PeptideMatch> getSeqMatch() { return mSeqMatch; }

  public ProteinID updateBestQVal(Double s) { if (mBestQVal==null || (s!=null && mBestQVal>s)) mBestQVal=s; return this; }

  public ProteinID put(String seq, PeptideMatch match)
  {
    if (mSeqMatch==null) mSeqMatch =HashMultimap.create();
    mSeqMatch.put(seq, match);

    return this;
  }
  public ProteinID putAll(String seq, Collection<PeptideMatch> matches)
  {
    if (seq!=null && Tools.isSet(matches))
      for (PeptideMatch match :matches) put(seq, match);

    return this;
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
        Float abund = (float )protein.getSeqMatch().values().size();
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
