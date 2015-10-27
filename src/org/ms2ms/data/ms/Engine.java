package org.ms2ms.data.ms;

import org.ms2ms.math.Stats;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/** Information about the search engine
 *
 * By Wen Yu on 9/20/2015.
 */
public class Engine implements Comparable<Engine>
{
  private String mName,            // the name
                 mCanonicalScore,  // the main score used to rank the entries from the same engine
                 mProb;            // A calibrated probability score for cross engine comparison
  public static final Engine MASCOT    = new Engine("Mascot", "IonScore");
  public static final Engine SEQUEST   = new Engine("Sequest", "XCorr");
  public static final Engine MSGF      = new Engine("MSGF+", "MSGFScore", "QValue");
  public static final Engine ANDROMEDA = new Engine("Maxquant-Andromeda", "Score");
  public static final Engine CRUX      = new Engine("Crux-Tide", "xcorr score");
  public static final Engine AMANDA    = new Engine("MS Amanda", "AmandaScore");

  // round out the parameters for the known engines
  static
  {

  }

  public Engine(String n)                     { mName=n; }
  public Engine(String n, String canonical)   { mName=n; mCanonicalScore=canonical; }
  public Engine(String n, String c, String q) { mName=n; mCanonicalScore=c; mProb=q; }

  public String getName()           { return mName; }
  public String getCanonicalScore() { return mCanonicalScore; }
  public String getProbScore()      { return mProb; }

  @Override
  public int compareTo(Engine o) { return mName.compareTo(o.getName()); }

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
}
