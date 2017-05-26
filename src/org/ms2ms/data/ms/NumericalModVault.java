package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.proteomics.mol.modification.unimod.UnimodManager;
import org.expasy.mzjava.proteomics.mol.modification.unimod.UnimodMod;
import org.ms2ms.utils.Tools;

import java.util.Collection;
import java.util.Collections;
import java.util.SortedMap;
import java.util.TreeMap;

/** keeper of delta mass, frequency and likely assignments
 *
 * Created by yuw on 5/26/17.
 */
public class NumericalModVault
{
  private String mVaultFile;
  private SortedMap<Double, Double>[] sAAMods=null;

  public NumericalModVault() { init(); }
  public NumericalModVault(String vault)
  {
    mVaultFile=vault; init();
  }
  public void init()
  {
    
    // setup the mod bank for each residues
    UnimodManager unimod = UnimodManager.getInstance();
    for (UnimodMod mod : unimod.getModificationList())
      if (Tools.isSet(mod.getSites()))
        for (String site : mod.getSites()) putAAMod(site, mod.getMolecularMass(), 1);
  }
  private void putAAMod(String aa, double mod, double prob)
  {
    if (sAAMods==null) sAAMods = new TreeMap[256];

    char A='*';
    if (aa.length()==1) A=aa.charAt(0);
    else if ("C-term".equals(aa)) A='$';
    else if ("N-term".equals(aa)) A='^';
    else
      System.out.print(aa);

    if (sAAMods[A]==null) sAAMods[A] = new TreeMap<Double, Double>();
    sAAMods[A].put(mod, sAAMods[A].get(mod)==null?prob:(prob+sAAMods[A].get(mod)));
  }

  public double bestScore(double mh, double delta, Tolerance tol, Collection<Integer> AAs)
  {
    if (AAs.size()==0) return 0;

    double d = 0.5d*(tol.getMax(mh)-tol.getMin(mh)), score=0;
    for (Integer A : AAs)
    {
      SortedMap<Double, Double> slice = sAAMods[A].subMap(delta-d, delta+d);
      if (Tools.isSet(slice) && Collections.max(slice.values())>score) score = Collections.max(slice.values());
    }
    return score;
  }
}
