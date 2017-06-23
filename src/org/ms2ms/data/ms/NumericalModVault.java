package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.proteomics.mol.modification.unimod.UnimodManager;
import org.expasy.mzjava.proteomics.mol.modification.unimod.UnimodMod;
import org.ms2ms.utils.TabFile;
import org.ms2ms.utils.Tools;

import java.util.*;

/** keeper of delta mass, frequency and likely assignments
 *
 * Created by yuw on 5/26/17.
 */
public class NumericalModVault
{
  private String mVaultFile="/Users/yuw/Apps/pipeline/Chorus/data/KnownModMasses.csv";
  private SortedMap<Double, Double>[] sAAMods=null, mKnown;
  private SortedMap<Double, Double> mDeltaProb;

  public NumericalModVault() { init(); }
  public NumericalModVault(String vault)
  {
    mVaultFile=vault; init();
  }
  public void init()
  {
    try
    {
      mDeltaProb = new TreeMap<>();
      TabFile vault = new TabFile(mVaultFile, TabFile.comma);
      while (vault.hasNext())
      {
        mKnown = putAAMod(mKnown, vault.get("Sites"), vault.getDouble("DeltaM"), Math.log10(vault.getDouble("Freq.mil.msms")));
        mDeltaProb.put(vault.getDouble("DeltaM"), Math.log10(vault.getDouble("Freq.mil.msms")));
      }
    }
    catch (Exception e) { e.printStackTrace(); }

    Collection<Double> checked = new HashSet<>();

    // setup the mod bank for each residues
    UnimodManager unimod = UnimodManager.getInstance();
    for (UnimodMod mod : unimod.getModificationList())
    {
      // any finding from the known?
      Map<Double, Double> slice = mDeltaProb!=null?mDeltaProb.subMap(mod.getMolecularMass()-0.01, mod.getMolecularMass()+0.01):null;
      double prob = Tools.isSet(slice) ? Collections.max(slice.values()):1;

      if (slice!=null) checked.addAll(slice.keySet());

      if (Tools.isSet(mod.getSites()))
        for (String site : mod.getSites())
          sAAMods = putAAMod(sAAMods, site, mod.getMolecularMass(), prob);
    }
    if (mDeltaProb!=null) mDeltaProb.keySet().removeAll(checked);
    if (Tools.isSet(mDeltaProb))
      for (Double d : mDeltaProb.keySet())
        sAAMods = putAAMod(sAAMods, "*", d, mDeltaProb.get(d));
  }
  private SortedMap<Double, Double>[] putAAMod(SortedMap<Double, Double>[] AAMods, String aa, double mod, double prob)
  {
    if (AAMods==null) AAMods = new TreeMap[256];

    char A='*';
    if (aa.length()==1) A=aa.charAt(0);
    else if ("C-term".equals(aa)) A='$';
    else if ("N-term".equals(aa)) A='^';
    else
      System.out.print(aa);

    if (AAMods[A]==null) AAMods[A] = new TreeMap<Double, Double>();
    AAMods[A].put(mod, AAMods[A].get(mod)==null?prob:(prob+AAMods[A].get(mod)));

    return AAMods;
  }

  public double bestScore(double mh, double delta, Tolerance tol, Collection<Integer> AAs)
  {
    double d = 0.5d*(tol.getMax(mh)-tol.getMin(mh)), score=0;
    if (Tools.isSet(AAs))
      for (Integer A : AAs)
      {
        if (delta-d>delta+d)
          System.out.println();
        SortedMap<Double, Double> slice = sAAMods[A].subMap(delta-d, delta+d);
        if (Tools.isSet(slice) && Collections.max(slice.values())>score) score = Collections.max(slice.values());
      }
    else
      for (int i=0; i<sAAMods.length; i++)
      {
        SortedMap<Double, Double> slice = sAAMods[i]!=null?sAAMods[i].subMap(delta-d, delta+d):null;
        if (Tools.isSet(slice) && Collections.max(slice.values())>score) score = Collections.max(slice.values());
      }

    return score;
  }
}