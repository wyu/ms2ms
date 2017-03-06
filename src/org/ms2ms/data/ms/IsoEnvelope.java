package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.ms2ms.algo.Isotopes;
import org.ms2ms.algo.Isotopics;
import org.ms2ms.algo.Peaks;
import org.ms2ms.algo.Similarity;
import org.ms2ms.utils.Tools;

import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;

/**
 * User: wyu
 * Date: Apr 5, 2011
 */
public class IsoEnvelope extends Peak
{
  private double mScore = 0d, mChargeScore = 0d;
  private List<Peak> mIsotopes = new ArrayList<Peak>(), mPredicted = new ArrayList<Peak>();

  public IsoEnvelope()        { super(); }
  public IsoEnvelope(Peak s) { super(s); }
  public IsoEnvelope(List<Peak> s, int charge)
  {
    super(s.get(0).getMz(), 1, charge);
    mPredicted.addAll(s);
  }
  public IsoEnvelope(double c12, int charge, double minri, double ai)
  {
    super(c12, ai, charge);

    double               limit = 0;
    List<Peak>           tmp = new ArrayList<>();
    Map<Integer, Long> formula = Isotopes.newFormulaMapByAveragine(c12 * charge);

    mPredicted = new ArrayList<>();

    // initialize the result
    mPredicted.add(new Peak(0.0, 1.0));

    Isotopes.calculate(tmp, mPredicted, formula, limit, charge);

    // move the predictions to the c12
    double offset = c12 - mPredicted.get(0).getMz(), min_ai = mPredicted.get(0).getIntensity()*minri*0.01;
    // scale the ai to the 1st c12
    ai /= mPredicted.get(0).getIntensity();

    Iterator<Peak> itr = mPredicted.iterator();
    while (itr.hasNext())
    {
      Peak pk = itr.next();
      pk.setMzAndCharge(pk.getMz()+offset, pk.getChargeList());
      pk.setIntensity(pk.getIntensity()*ai);
      if (pk.getIntensity() < min_ai) itr.remove();
    }
  }
  public IsoEnvelope(double c12, int charge, double minri, double ai, Isotopics iso)
  {
    super(c12, ai, charge);

    double               limit = 0;
    List<Peak>             tmp = new ArrayList<>();
    Map<Integer, Long> formula = iso.newFormulaMapByAveragine(c12 * charge);

    mPredicted = new ArrayList<>();

    // initialize the result
    mPredicted.add(new Peak(0.0, 1.0));

    iso.calculate(tmp, mPredicted, formula, limit, charge);

    // move the predictions to the c12
    double offset = c12 - mPredicted.get(0).getMz(), min_ai = mPredicted.get(0).getIntensity()*minri*0.01;
    // scale the ai to the 1st c12
    ai /= mPredicted.get(0).getIntensity();

    Iterator<Peak> itr = mPredicted.iterator();
    while (itr.hasNext())
    {
      Peak pk = itr.next();
      pk.setMzAndCharge(pk.getMz()+offset, pk.getChargeList());
      pk.setIntensity(pk.getIntensity()*ai);
      if (pk.getIntensity() < min_ai) itr.remove();
    }
  }

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
    return Similarity.dp(A, B)>0.9;
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
