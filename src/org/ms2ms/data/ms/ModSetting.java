package org.ms2ms.data.ms;

import org.ms2ms.Disposable;
import org.ms2ms.algo.Peaks;
import org.ms2ms.algo.Peptides;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Tools;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 // http://www.matrixscience.com/help/search_field_help.html
 // http://www.matrixscience.com/blog/ion-series-for-ethcd.html

 * Created by yuw on 8/20/16.
 */
public class ModSetting implements Disposable
{
  public enum eFrag { hcd, cid, ethcd }

  public double[]          masses;
  public char[]            charges, sumZ, sumH2O, sumNH3;
  public List<ModLocation> mods;

  @Deprecated
  public ModSetting(double[] m, List<ModLocation> mod)           { masses=m; mods=mod; }
  public ModSetting(int size, int loc, double increment)         { init(size,loc,increment); }
  public ModSetting(double[] m, char[] z, List<ModLocation> mod) { masses=m; charges=z; mods=mod; }

  // produce the mass ladder considering the AA residue, existing mods and incremental mass at loc
  public ModSetting(char[] peptide, Map<Integer, Double> Mods, int loc, double increment, ResidueBase base)
  {
    // other ion series? 20180317.
    if (peptide!=null && peptide.length>0)
    {
      Float incre_z = 0f;
      init(peptide.length, loc, increment);

      if (increment!=0 && base.getModCharge()!=null)
      {
        Map<Double, Float> slice = base.getModCharge().subMap(increment-0.01, increment+0.01);
        if (Tools.isSet(slice)) incre_z = Stats.sumFloats(slice.values())*10f;
      }
      char h2o=0, nh3=0, z=0;
      for (int i=0; i<peptide.length; i++)
      {
        Double dM = (Mods!=null&&Mods.containsKey(i)?Mods.get(i):0d);
        h2o += "STED".indexOf(peptide[i]>=0?1:0); nh3 += "RKNQ".indexOf(peptide[i]>=0?1:0);
        z += "RK".indexOf(peptide[i]>=0?10:0); z += "H".indexOf(peptide[i]>=0?5:0);
        masses[ i] = base.getAAs()[peptide[i]] + (i==loc?increment:0d) + dM;
        sumH2O[ i] = h2o; sumNH3[i] = nh3; sumZ[i]=z;
        charges[i] = (char )Math.round(("KR".indexOf(peptide[i])>=0?10:("H".indexOf(peptide[i])>=0?5:0)) +
            (i==loc?incre_z:0) + (base.getModCharge().containsKey(dM)?base.getModCharge().get(dM)*10f:0f));
      }
    }
//    Float incre_z = 0f;
//    ModSetting mm = new ModSetting(peptide.length, loc, increment);
//
//    if (increment!=0 && getBase().getModCharge()!=null)
//    {
//      Map<Double, Float> slice = getBase().getModCharge().subMap(increment-0.01, increment+0.01);
//      if (Tools.isSet(slice)) incre_z = Stats.sumFloats(slice.values())*10f;
//    }
//    for (int i=0; i<peptide.length; i++)
//    {
//      Double dM = (mods!=null&&mods.containsKey(i)?mods.get(i):0d);
//      mm.masses[ i] = AAs[peptide[i]] + (i==loc?increment:0d) + dM;
//      mm.charges[i] = (char )Math.round(("KR".indexOf(peptide[i])>=0?10:("H".indexOf(peptide[i])>=0?5:0)) +
//          (i==loc?incre_z:0) + (getBase().getModCharge().containsKey(dM)?getBase().getModCharge().get(dM)*10f:0f));
//    }
  }
  // calc the y and other c-term series from a prev stop
  public double[] calcYs(int i, boolean multiZ, double prev)
  {
    int  x=masses.length-1, k=x-i;
    prev +=masses[k];
    // add other series if necessary
    double[] mm = new double[] {
        prev,                                      // simple y-series
        sumH2O[x]-sumH2O[k]>0?prev-Peptides.H2O:0, // loss of 18 (H2O) if STED in the fragment
        sumNH3[x]-sumNH3[k]>0?prev-Peptides.NH3:0  // loss of 17 (NH3) if RKNQ in the fragment
    };
    return multiZ&&sumZ[x]>sumZ[k]?Peaks.MH2Mzs(1,Math.round(0.1f*(sumZ[x]-sumZ[k]))):mm;
  }
  public double[] calcBs(int i, boolean multiZ, double prev)
  {
    prev +=masses[i];
    // add other series if necessary
    double[] mm = new double[] {
        prev,                                      // simple b-series
        prev-Peptides.CO,                          // a-series
        sumH2O[i]>0?prev-Peptides.H2O:0,           // loss of 18 (H2O) if STED in the fragment
        sumNH3[i]>0?prev-Peptides.NH3:0            // loss of 17 (NH3) if RKNQ in the fragment
    };
    return multiZ&&sumZ[i]>0?Peaks.MH2Mzs(1,Math.round(0.1f*(sumZ[i]))):mm;
  }
//  public double[] calcInternals(int i, boolean multiZ)
//  {
//    prev +=masses[i];
//    // add other series if necessary
//    double[] mm = new double[] {
//        prev,                                      // simple b-series
//        prev-Peptides.CO,                          // a-series
//        sumH2O[i]>0?prev-Peptides.H2O:0,           // loss of 18 (H2O) if STED in the fragment
//        sumNH3[i]>0?prev-Peptides.NH3:0            // loss of 17 (NH3) if RKNQ in the fragment
//    };
//    return multiZ&&sumZ[i]>0?Peaks.MH2Mzs(1,Math.round(0.1f*(sumZ[i]))):mm;
//  }
  private void init(int size, int loc, double increment)
  {
    masses = new double[size];
    charges= new char[  size];
    sumZ   = new char[  size];
    sumH2O = new char[  size];
    sumNH3 = new char[  size];

    mods = new ArrayList<>(); mods.add(new ModLocation(loc, increment));
  }

  @Override
  public void dispose()
  {
    masses = null;
    charges= null;
    sumZ   = null;
    sumH2O = null;
    sumNH3 = null;

    mods = (List )Tools.dispose(mods);
  }
}
