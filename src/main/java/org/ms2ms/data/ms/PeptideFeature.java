package org.ms2ms.data.ms;

import com.google.common.base.Optional;
import org.expasy.mzjava.core.mol.Mass;
import org.expasy.mzjava.proteomics.mol.AminoAcid;
import org.expasy.mzjava.proteomics.mol.modification.ModAttachment;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.ms.ident.ModificationMatch;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.expasy.mzjava.proteomics.ms.ident.PeptideProteinMatch;
import org.ms2ms.algo.Isobarics;
import org.ms2ms.data.Binary;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

/**
 * ** Copyright 2014-2015 ms2ms.org
 * <p/>
 * Description: A peptide ID with quantitative info
 * <p/>
 * Author: wyu
 * Date:   3/26/15
 */
public class PeptideFeature extends PeptideMatch implements Binary
{
  public static final String SCORE = "Score";

  private ProteinID mProteinID;
  private Collection<PeptideMatch> mMatches;
  private Map<String, Double> mExptAbundance = new HashMap<>();;
  private int mCharge;
  private String mTitle;
  private Double mRT, mMZ;

  public PeptideFeature() { super(); }
  public PeptideFeature(String s) { super(s); }
  public PeptideFeature(PeptideMatch s)
  {
    super(s);
  }

  public ProteinID          getProteinID()         { return mProteinID; }
  public String             getTitle()             { return mTitle; }
  public int                getCharge()            { return mCharge; }
  public Double             getRT()                { return mRT; }
  public Double             getMz()                { return mMZ; }
  public Collection<PeptideMatch> getMatches()     { return mMatches; }
  public Double             getAbundance(String s) { return mExptAbundance.get(s); }
  public Collection<String> getExperiments()       { return mExptAbundance!=null?mExptAbundance.keySet():null; }

  public PeptideFeature setProteinID(ProteinID        s){ mProteinID=s; return this; }
  public PeptideFeature setAbundance(String t, Double s){ mExptAbundance.put(t, s); return this; }
  public PeptideFeature setCharge(int s)                { mCharge=s; return this; }
  public PeptideFeature setRT(Double s)                 { mRT=s; return this; }
  public PeptideFeature setMz(Double s)                 { mMZ=s; return this; }
  public PeptideFeature setTitle(String s)              { mTitle=s; return this; }
  public PeptideFeature addAbundance(String t, Double s)
  {
    mExptAbundance.put(t, getAbundance(t)!=null?s+getAbundance(t):s);
    return this;
  }
  public PeptideFeature addMatch(PeptideMatch s)
  {
    if (mMatches==null) mMatches = new ArrayList<>();
    mMatches.add(s);
    return this;
  }

  @Override
  public int hashCode()
  {
    return Tools.hashCode(getTitle(), getProteinID().getAccession());
  }
  @Override
  public String toString()
  {
    return super.toString();
  }

  @Override
  public int compareTo(PeptideMatch o)
  {
    int OK = super.compareTo(o);
    if (OK==0 && hasScore(SCORE) && o.hasScore(SCORE))
        OK = Double.compare(getScore(SCORE), o.getScore(SCORE));
    if (OK==0 && Tools.isSet(mExptAbundance) && Tools.isSet(((PeptideFeature) o).mExptAbundance))
        OK = Double.compare(Stats.sum(mExptAbundance.values()), Stats.sum(((PeptideFeature) o).mExptAbundance.values()));
    return OK;
  }

  /** create a new peptide feature from the output of 'multiplierz'
   *
   *
   * @param title is the title of the experiment
   * @param set containing a particular row of the scan-level summary of a TMT experiment
   * @param expts are the cols for the sample channels
   * @param base is the abundance of the control experiment. Not check against 0
   * @return a new peptide feature
   */
  public static PeptideFeature fromMultiplierZ(String title, Map<String, String> set, Map<String, String> expts, double base)
  {
    PeptideFeature peptide = new PeptideFeature(set.get("Peptide Sequence"));
    peptide.setCharge(new Integer(set.get("Charge"))).setTitle(title);
    peptide.setMassDiff(Stats.toDouble(set.get("Delta")));
    peptide.addScore(PeptideFeature.SCORE, Stats.toDouble(set.get("Peptide Score")));
    peptide.addProteinMatch(new PeptideProteinMatch(set.get("Accession Number"), Optional.of("null"),
            Optional.of(set.get("Preceding Residue")), Optional.of(set.get("Following Residue")), PeptideProteinMatch.HitType.UNKNOWN));
    //S15: Phospho; S22: Phospho
    String[] mods = Strs.split(set.get("Variable Modifications"), ';', true);
    if (Tools.isSet(mods))
      for (String mod : mods)
      {
        String[] parsed = Strs.split(mod, ':', true);
        if (parsed!=null && parsed.length>1 && parsed[0].length()>1 && Stats.toInt(parsed[0].substring(1))!=null)
          peptide.addModificationMatch(Stats.toInt(parsed[0].substring(1))-1,
              new ModificationMatch(new Modification(parsed[1]+"@"+parsed[0].charAt(0), Mass.ZERO),
                  AminoAcid.valueOf(parsed[0]), Integer.valueOf(parsed[1]), ModAttachment.SIDE_CHAIN)); // default to side-chain for now
      }
    if (Tools.isSet(expts))
      for (String expt : expts.keySet())
        peptide.setAbundance(expts.get(expt), Stats.toDouble(set.get(expt))!=null?Stats.toDouble(set.get(expt))/base:0d);

    // set the basis
    return peptide.setAbundance(Isobarics.BASIS, base);
  }

  @Override
  public void write(DataOutput ds) throws IOException
  {
//    private ProteinID mProteinID;
//    private Collection<PeptideMatch> mMatches;
//    private Map<String, Double> mExptAbundance = new HashMap<>();;
//    private int mCharge;
//    private String mTitle;
//    private Double mRT, mMZ;
  }

  @Override
  public void read(DataInput ds) throws IOException
  {

  }
}
