package org.ms2ms.data.ms;

import gnu.trove.map.TObjectDoubleMap;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import javax.annotation.Nonnull;
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
                 mProb;            // A calibrated probability score for cross engine comparison
  public static final Engine MASCOT    = new Engine("Mascot", "IonScore");
  public static final Engine SEQUEST   = new Engine("Sequest", "XCorr");
  public static final Engine MSGF      = new Engine("MSGF+", "MSGFScore", "QValue");
  public static final Engine ANDROMEDA = new Engine("Maxquant-Andromeda", "Score");
  public static final Engine CRUX      = new Engine("Crux-Tide", "xcorr score");
  public static final Engine AMANDA    = new Engine("MS Amanda", "AmandaScore");
  public static final Engine XTANDEM   = new Engine("X!Tandem", "X\\!Tandem:hyperscore", "X\\!Tandem:expect");
  public static final Engine OMSSA     = new Engine("OMSSA", "X\\!Tandem:hyperscore", "X\\!Tandem:expect");

  // round out the parameters for the known engines
  static
  {

  }

  public Engine(String n)                     { mName=n; }
  public Engine(String n, String canonical)   { mName=n; mCanonicalScore=canonical; }
  public Engine(String n, String canonical, boolean desending)   { mName=n; mCanonicalScore=canonical; mIsDesending=desending; }
  public Engine(String n, String c, String q) { mName=n; mCanonicalScore=c; mProb=q; }
  public Engine(String n, String c, String q, boolean desending) { mName=n; mCanonicalScore=c; mProb=q; mIsDesending=desending; }

  public String getName()           { return mName; }
  public String getCanonicalScore() { return mCanonicalScore; }
  public String getProbScore()      { return mProb; }
  public boolean isDesending()      { return mIsDesending; }

  public Engine isDesending(boolean s) { mIsDesending=s; return this; }

  @Override
  public int compareTo(Engine o) { return mName.compareTo(o.getName()); }

  public int compareTo(@Nonnull TObjectDoubleMap<String> A, @Nonnull TObjectDoubleMap<String> B)
  {
    return (isDesending()?
      Stats.compareTo(B.get(getCanonicalScore()), A.get(getCanonicalScore())):
      Stats.compareTo(A.get(getCanonicalScore()), B.get(getCanonicalScore())));
  }
  //** utilitiy functions **//
  public static Engine valueOf(String s)
  {
    if      (Strs.equalsIgnoreCase(s,    MASCOT.getName())) return MASCOT;
    else if (Strs.equalsIgnoreCase(s,   SEQUEST.getName())) return SEQUEST;
    else if (Strs.equalsIgnoreCase(s,      MSGF.getName())) return MSGF;
    else if (Strs.equalsIgnoreCase(s, ANDROMEDA.getName())) return ANDROMEDA;

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
