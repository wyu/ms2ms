package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.ms2ms.algo.Peaks;
import org.ms2ms.algo.Similarity;
import org.ms2ms.utils.Tools;

import java.util.List;
import java.util.ArrayList;

/**
 * User: wyu
 * Date: Apr 5, 2011
 */
public class IsoEnvelope
{
  private double mz = 0.0D;
  private double intensity = 0.0D;
  private int charge = 0;
//  private double mass = 0.0D;

  private double mScore = 0d, mChargeScore = 0d;
  private List<Peak> mIsotopes = new ArrayList<Peak>(), mPredicted = new ArrayList<Peak>();

  public IsoEnvelope()        { super(); }
//  public IsoEnvelope(Peak s)
//  {
//    init(mz, s.getIntensity(), s.getCharge());
//  }
  public IsoEnvelope(List<Peak> s, int z)
  {
    init(s.get(0).getMz(),1,z);
    // use the list directly since we are disposing it anywhere. WYU 20180118
//    mPredicted=s;
    mPredicted.addAll(s);
  }
  public IsoEnvelope(double c12, double ai, int z) { init(c12, ai, z); }
//  public IsoEnvelope(double c12, int z, double minri, double ai)
//  {
//    init(c12, ai, z);
//
//    double               limit = 0;
//    List<Peak>           tmp = new ArrayList<>();
//    Map<Integer, Long> formula = Isotopes.newFormulaMapByAveragine(c12 * charge);
//
//    mPredicted = new ArrayList<>();
//
//    // initialize the result
//    mPredicted.add(new Peak(0.0, 1.0));
//
//    Isotopes.calculate(tmp, mPredicted, formula, limit, charge);
//
//    // move the predictions to the c12
//    double offset = c12 - mPredicted.get(0).getMz(), min_ai = mPredicted.get(0).getIntensity()*minri*0.01;
//    // scale the ai to the 1st c12
//    ai /= mPredicted.get(0).getIntensity();
//
//    Iterator<Peak> itr = mPredicted.iterator();
//    while (itr.hasNext())
//    {
//      Peak pk = itr.next();
//      pk.setMzAndCharge(pk.getMz()+offset, pk.getChargeList());
//      pk.setIntensity(pk.getIntensity()*ai);
//      if (pk.getIntensity() < min_ai) itr.remove();
//    }
//  }
//  public IsoEnvelope(double c12, int z, double minri, double ai, Isotopics iso)
//  {
//    init(c12, ai, z);
//
//    double               limit = 0;
//    List<Peak>             tmp = new ArrayList<>();
//    Map<Integer, Long> formula = iso.newFormulaMapByAveragine(c12 * charge);
//
//    mPredicted = new ArrayList<>();
//
//    // initialize the result
//    mPredicted.add(new Peak(0.0, 1.0));
//
//    iso.calculate(tmp, mPredicted, formula, limit, charge);
//
//    // move the predictions to the c12
//    double offset = c12 - mPredicted.get(0).getMz(), min_ai = mPredicted.get(0).getIntensity()*minri*0.01;
//    // scale the ai to the 1st c12
//    ai /= mPredicted.get(0).getIntensity();
//
//    Iterator<Peak> itr = mPredicted.iterator();
//    while (itr.hasNext())
//    {
//      Peak pk = itr.next();
//      pk.setMzAndCharge(pk.getMz()+offset, pk.getChargeList());
//      pk.setIntensity(pk.getIntensity()*ai);
//      if (pk.getIntensity() < min_ai) itr.remove();
//    }
//  }
  public IsoEnvelope init(double m, double ai, int z)
  {
    mz=m; intensity=ai; charge=z;
    return this;
  }
  public double getMz() { return mz; }
  public double getIntensity() { return intensity; }
  public int getCharge() { return charge; }
  public IsoEnvelope setMz(double s) { mz=s; return this; }
  public IsoEnvelope setIntensity(double s) { intensity=s; return this; }
  public IsoEnvelope setCharge(int s) { charge=s; return this; }

  public IsoEnvelope setChargeScore(double s) { mChargeScore = s; return this; }
  public IsoEnvelope setScore(double s) { mScore = s; return this; }
  public IsoEnvelope addIsotope(Peak ion)
  {
    if (!mIsotopes.contains(ion)) mIsotopes.add(ion); return this;
  }
  public IsoEnvelope addPredicted(Peak ion)
  {
    mPredicted.add(ion); return this;
  }
  public boolean contains(IsoEnvelope other, Tolerance tol)
  {
    if (other == null) return false;

    List<Peak> A = new ArrayList<Peak>(), B = new ArrayList<Peak>();
    for (Peak k : other.mIsotopes)
    {
      boolean found = false;
      for (Peak i : mIsotopes)
        if (tol.withinTolerance(i.getMz(), k.getMz()))
        { found = true; A.add(k); B.add(i); break; }
      //if (!found) return false;
      if (!found) { A.add(k); B.add(new Peak(k.getMz(), 0d, 1)); }
    }
    return Similarity.dp(A, B, false)>0.9;
  }
  public double      getScore()          { return mScore; }
  public double      getChargeScore()    { return mChargeScore; }
  public List<Peak>  getIsotopes()       { return mIsotopes; }
  public List<Peak>  getPredicted()      { return mPredicted; }
  public      Peak   get(int s)          { return s >= 0 && s < mIsotopes.size() ? mIsotopes.get(s) : null; }
  public      Peak   getPredicted(int s) { return s >= 0 && s < mPredicted.size() ? mPredicted.get(s) : null; }

  public String toString()
  {
    double max = Peaks.getBasePeak(getIsotopes()).getIntensity();
    StringBuffer buf = new StringBuffer();

    buf.append("m/z" + Tools.d2s(getMz(), 4) + "@+" + getCharge() + "&" + Tools.d2s(getIntensity(), 0) + ", ");
    buf.append("zscore/score: " + Tools.d2s(getChargeScore(), 1) + "/" + Tools.d2s(getScore(), 2) + ">>");
    for (Peak pt : getIsotopes())
      buf.append(Tools.d2s(pt.getMz(), 3) + "," + Tools.d2s(pt.getIntensity() * 100d / max, 1) + "; ");

    return buf.toString();
  }
}
