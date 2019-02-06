package org.ms2ms.data.ms;

import com.google.common.collect.Multimap;
import gnu.trove.map.TObjectDoubleMap;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import org.ms2ms.algo.PSMs;
import org.ms2ms.io.PsmReaders;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/** Information about the search engine
 *
 * By Wen Yu on 9/20/2015.
 */
public class Engine implements Comparable<Engine>
{
  private boolean mIsDesending = true;
  private String mName=null,            // the name
                 mPsiName=null,
                 mVersion="0.0",
                 mCanonicalScore=null,  // the main score used to rank the entries from the same engine
                 mDeltaScore=PSMs.SCR_DELTA,
                 mProb=null,            // A calibrated probability score for cross engine comparison
                 mPsiCV=null;
  private String[] mPercolators = null;

  public static final Engine AMANDA    = new Engine("MS Amanda", "AmandaScore").setPercolators(
      "SpecId", "Label","ScanNr","ExpMass","CalcMass","Charge1", "Charge2", "Charge3", "Charge4", "Mass", "PepLen", "enzC", "dM", "absdM",
      "AmandaScore","WeightedProbability","yb", "Peptide", "Proteins");
  public static final Engine ANDROMEDA = new Engine("Maxquant-Andromeda", "Score");
  public static final Engine BYONIC    = new Engine("Byonic", "Score");
  public static final Engine CRUX      = new Engine("Crux-Tide", "xcorr_score").setPercolators(
      "SpecId", "Label", "ScanNr","ExpMass","CalcMass", "Charge1", "Charge2", "Charge3", "Charge4", "Mass", "PepLen", "enzN", "enzC", "dM", "absdM",
      PSMs.SCR_DELTA,"sp_score","xcorr_score","delta_cn","IonFrac","lnrSp", "yb", "Peptide", "Proteins");
  public static final Engine MASCOT    = new Engine("Mascot", "IonScore");
  public static final Engine MAXQUANT    = new Engine("MaxQuant", "Score");
  public static final Engine MSGF      = new Engine("MSGF+", "MSGFScore", "QValue").setPercolators(
    "SpecId","Label","ScanNr","ExpMass","CalcMass","Charge1","Charge2","Charge3","Charge4","Mass","PepLen","enzN","enzC","dM","absdM",
    "DeNovoScore","IsotopeError","MSGFScore","SpecEValue","EValue","QValue","PepQValue", "yb","Peptide","Proteins");
  public static final Engine MYRIMATCH = new Engine("MyriMatch", "MyriMatch:MVH").setPercolators(
      "SpecId","Label","ScanNr","ExpMass","CalcMass","Charge1","Charge2","Charge3","Charge4","Mass","PepLen","enzN","enzC","dM","absdM",
      "xcorr","MyriMatch:mzFidelity","MyriMatch:MVH", "yb","Peptide","Proteins");
  public static final Engine OMSSA     = new Engine("OMSSA", "dbEVal").setPercolators(
      "SpecId","Label","ScanNr","ExpMass","CalcMass","Charge1","Charge2","Charge3","Charge4","Mass","PepLen","enzN","enzC","dM","absdM",
      "OMSSA:EVal","OMSSA:PVal","dbEVal","dbPVal","OMSSA:NIST", "yb","Peptide","Proteins");
  public static final Engine PEAKS     = new Engine("PEAKS", "PeakScore");
  public static final Engine PERCOLATOR= new Engine("Percolator", "qVal", "PEP");
  public static final Engine PPILOT    = new Engine("ProteinPilot", "Sc", "Conf");
  public static final Engine SEQUEST   = new Engine("Sequest", "XCorr");
  public static final Engine XTANDEM   = new Engine("X!Tandem", "X\\!Tandem:hyperscore").setPercolators(
      "SpecId","Label","ScanNr","ExpMass","CalcMass","Charge1","Charge2","Charge3","Charge4","Mass","PepLen","enzN","enzC","dM","absdM",
      "X\\!Tandem:hyperscore", "X\\!Tandem:expect", "yb","Peptide","Proteins").setPsiName("X\\!Tandem").setPsiCV("MS:1001476");

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
  public String getPsiName()        { return mPsiName; }
  public String getCanonicalScore() { return mCanonicalScore; }
  public String getProbScore()      { return mProb; }
  public String getVersion()        { return mVersion; }
  public String getPsiCV()          { return mPsiCV; }
  public boolean isDesending()      { return mIsDesending; }
  public String[] getPercolators()  { return mPercolators; }

  public Engine setDeltaScore(    String s) { mDeltaScore=s;  return this; }
  public Engine setPsiCV(         String s) { mPsiCV=s;       return this; }
  public Engine setPsiName(       String s) { mPsiName=s;     return this; }
  public Engine isDesending(     boolean s) { mIsDesending=s; return this; }
  public Engine setPercolators(String... s) { mPercolators=s; return this; }

  @Override
  public int compareTo(Engine o) { return mName.compareTo(o.getName()); }

  public int compareTo(TObjectDoubleMap<String> A, TObjectDoubleMap<String> B)
  {
    return (isDesending()?
      Stats.compareTo(B.get(getCanonicalScore()), A.get(getCanonicalScore())):
      Stats.compareTo(A.get(getCanonicalScore()), B.get(getCanonicalScore())));
  }
  public Multimap<SpectrumIdentifier, PeptideMatch> readPSMs(String s, int lowest_rank) throws IOException
  {
    if      (Strs.equals(     MSGF.getName(), mName)) return PsmReaders.readMSGFplus(s, '\t', lowest_rank);
    else if (Strs.equals(    OMSSA.getName(), mName)) return PsmReaders.readOMSSA(s);
    else if (Strs.equals(     CRUX.getName(), mName)) return PsmReaders.readCrux(s, lowest_rank);
    else if (Strs.equals(MYRIMATCH.getName(), mName)) return PsmReaders.readMyriMatch(s, lowest_rank);
    else if (Strs.equals(  XTANDEM.getName(), mName)) return PsmReaders.readXTandem(s, lowest_rank);

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

  public static Engine[] values()
  {
    return new Engine[] {AMANDA,ANDROMEDA,BYONIC,CRUX,MASCOT,MAXQUANT,MSGF,MYRIMATCH,OMSSA,PEAKS,PERCOLATOR,PPILOT,SEQUEST,XTANDEM};
  }
  public static String[] names()
  {
    String[] names = new String[values().length];
    for (int i=0; i< values().length; i++) names[i]=values()[i].getName();

    return names;
  }

  @Override
  public String toString()
  {
    return getName() + ": " + getCanonicalScore() + "/" + getProbScore();
  }
}
