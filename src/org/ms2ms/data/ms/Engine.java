package org.ms2ms.data.ms;

import com.google.common.collect.Multimap;
import gnu.trove.map.TObjectDoubleMap;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import org.ms2ms.algo.PSMs;
import org.ms2ms.io.PsmReaders;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.TabFile;
import org.ms2ms.utils.Tools;

import javax.annotation.Nonnull;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/** Information about the search engine
 *
 * By Wen Yu on 9/20/2015.
 */
public class Engine implements Comparable<Engine>
{
  private boolean mIsDesending = true;
  private String mName,            // the name
                 mCanonicalScore,  // the main score used to rank the entries from the same engine
                 mDeltaScore=PSMs.SCR_DELTA,
                 mProb;            // A calibrated probability score for cross engine comparison
  private String[] mPercolators = null;

  public static final Engine AMANDA    = new Engine("MS Amanda", "AmandaScore").setPercolators(
      "SpecId", "Label", "ScanNr", "Charge1", "Charge2", "Charge3", "Charge4", "Mass", "PepLen", "enzC", "dM", "absdM",
      "AmandaScore","WeightedProbability","yb", "Peptide", "Proteins");
  public static final Engine ANDROMEDA = new Engine("Maxquant-Andromeda", "Score");
  public static final Engine BYONIC    = new Engine("Byonic", "");
  public static final Engine CRUX      = new Engine("Crux-Tide", "xcorr score").setPercolators(
      "SpecId", "Label", "ScanNr", "Charge1", "Charge2", "Charge3", "Charge4", "Mass", "PepLen", "enzN", "enzC", "dM", "absdM",
      PSMs.SCR_DELTA,"sp_score","xcorr_score","IonFrac","lnrSp", "yb", "Peptide", "Proteins");
  public static final Engine MASCOT    = new Engine("Mascot", "IonScore");
  public static final Engine MSGF      = new Engine("MSGF+", "MSGFScore", "QValue").setPercolators(
    "SpecId","Label","ScanNr","Charge1","Charge2","Charge3","Charge4","Mass","PepLen","enzN","enzC","dM","absdM",
    "DeNovoScore","IsotopeError","MSGFScore","SpecEValue","EValue","QValue","PepQValue", "yb","Peptide","Proteins");
  public static final Engine MYRIMATCH = new Engine("MyriMatch", "", "").setPercolators();
  public static final Engine OMSSA     = new Engine("OMSSA", "X\\!Tandem:hyperscore", "X\\!Tandem:expect").setPercolators();
  public static final Engine PEAKS     = new Engine("PEAKS", "PeakScore");
  public static final Engine PPILOT    = new Engine("ProteinPilot", "Sc", "Conf");
  public static final Engine SEQUEST   = new Engine("Sequest", "XCorr");
  public static final Engine XTANDEM   = new Engine("X!Tandem", "X\\!Tandem:hyperscore", "X\\!Tandem:expect");

  // round out the parameters for the known engines
  static
  {

  }

  public Engine(String n)                     { mName=n; }
  public Engine(String n, String canonical)   { mName=n; mCanonicalScore=canonical; }
  public Engine(String n, String canonical, boolean desending)   { mName=n; mCanonicalScore=canonical; mIsDesending=desending; }
  public Engine(String n, String c, String q) { mName=n; mCanonicalScore=c; mProb=q; }
  public Engine(String n, String c, String q, boolean desending) { mName=n; mCanonicalScore=c; mProb=q; mIsDesending=desending; }

  public String getDeltaScore()     { return mDeltaScore; }
  public String getName()           { return mName; }
  public String getCanonicalScore() { return mCanonicalScore; }
  public String getProbScore()      { return mProb; }
  public boolean isDesending()      { return mIsDesending; }
  public String[] getPercolators()  { return mPercolators; }

  public Engine setDeltaScore(    String s) { mDeltaScore=s; return this; }
  public Engine isDesending(     boolean s) { mIsDesending=s; return this; }
  public Engine setPercolators(String... s) { mPercolators=s; return this; }

  @Override
  public int compareTo(Engine o) { return mName.compareTo(o.getName()); }

  public int compareTo(@Nonnull TObjectDoubleMap<String> A, @Nonnull TObjectDoubleMap<String> B)
  {
    return (isDesending()?
      Stats.compareTo(B.get(getCanonicalScore()), A.get(getCanonicalScore())):
      Stats.compareTo(A.get(getCanonicalScore()), B.get(getCanonicalScore())));
  }
  public Multimap<SpectrumIdentifier, PeptideMatch> readPSMs(String s) throws IOException
  {
    if      (Strs.equals( MSGF.getName(), mName)) return PsmReaders.readMSGFplus(s, '\t');
    else if (Strs.equals(OMSSA.getName(), mName)) return PsmReaders.readOMSSA(s);

    return null;
  }
  //** utilitiy functions **//
  public static Engine valueOf(String s)
  {
    if      (Strs.equalsIgnoreCase(s, ANDROMEDA.getName())) return ANDROMEDA;
    else if (Strs.equalsIgnoreCase(s,    MASCOT.getName())) return MASCOT;
    else if (Strs.equalsIgnoreCase(s,   SEQUEST.getName())) return SEQUEST;
    else if (Strs.equalsIgnoreCase(s,      MSGF.getName())) return MSGF;
    else if (Strs.equalsIgnoreCase(s, MYRIMATCH.getName())) return MYRIMATCH;
    else if (Strs.equalsIgnoreCase(s,   XTANDEM.getName())) return XTANDEM;
    else if (Strs.equalsIgnoreCase(s,     OMSSA.getName())) return OMSSA;
    else if (Strs.equalsIgnoreCase(s,      CRUX.getName())) return CRUX;
    else if (Strs.equalsIgnoreCase(s,    AMANDA.getName())) return AMANDA;

    return null;
  }

  public static Engine[] valuesOf(String... s)
  {
    if (!Tools.isSet(s)) return null;

    List<Engine> engines = new ArrayList<>();
    for (String name : s)
      Tools.add(engines, valueOf(name));

    return engines.toArray(new Engine[] {});
  }

  @Override
  public String toString()
  {
    return getName() + ": " + getCanonicalScore() + "/" + getProbScore();
  }
}
