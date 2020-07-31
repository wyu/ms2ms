package org.ms2ms.data.ms;

import org.ms2ms.data.Point;
import org.ms2ms.math.Points;
import org.ms2ms.utils.Tools;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;

/** A representation of LcMs feature from MaxQuant or similar program.
 *
 */
public class LcMsFeature extends LcMsPoint
{
  public enum criteria { snr, dp, ms1, dpms1, inbound, suburb, outer, top }

  private boolean mInMS1=false;
  private int    mApexPos, mPointWidth, mPresentN;
  private double mArea, mApexRT, mLower, mUpper;
  private double mAbundance;
  private double mApex, mDeriv1st, mSNR;
  private double mMzStdev=Double.NaN, mMzStdevEx=Double.NaN;
  private double mInitialCentroidRt, mSimilarity, mExclusivity;
  private criteria mPeakSelection;

  public static class AreaDesendComparator implements Comparator<LcMsFeature> {
    public int compare(LcMsFeature o1, LcMsFeature o2) {
      return o1 != null && o2 != null ? Double.compare(o2.getArea(), o1.getArea()) : 0;
    }
  }
  public static class SimilarityDesendComparator implements Comparator<LcMsFeature> {
    public int compare(LcMsFeature o1, LcMsFeature o2) {
      return o1 != null && o2 != null ? Double.compare(o2.getSimilarity(), o1.getSimilarity()) : 0;
    }
  }

  public LcMsFeature() {};
  public LcMsFeature(double mz, double rt)
  {
  }
  public LcMsFeature(double rt, double ai, double d1)
  {
    setX(rt); setY(ai); mDeriv1st=d1;
  }
  public LcMsFeature(List<LcMsPoint> pts)
  {
    super();
    if (Tools.isSet(pts))
    {
      setX(Points.centroid(pts));
      if (Points.basePoint(pts)!=null) setY(Points.basePoint(pts).getY());
      setArea(Points.sumY(pts));
    }
  }

  public LcMsFeature(Point s)
  {
    if (s!=null) { setX(s.getX()); setY(s.getY()); }
  }

  public double getAbundance() { return mAbundance; }
  public double getArea()      { return mArea; }
  public double getApex()      { return mApex; }
  public double getApexRT()    { return mApexRT; }
  public int    getApexPos()   { return mApexPos; }
  public int    getPointWidth() { return mPointWidth; }
  public int    isPresentIn()   { return mPresentN; }
  public double getSNR()        { return mSNR; }
  public double getSimilarity() { return mSimilarity; }
  public double getExclusivity() { return mExclusivity; }
  public double getLower() { return mLower; }
  public double getUpper() { return mUpper; }
  public boolean inMS1() { return mInMS1; }

  public double getMzStdev()   { return mMzStdev; }
  public double getMzStdevEx() { return mMzStdevEx; }
  public double getInitialCentroidRt()      { return mInitialCentroidRt; }
  public criteria wasBasedOn() { return mPeakSelection; }

  public LcMsFeature isPresentIn(int s) { mPresentN=s; return this; }
  public LcMsFeature setApexPos(int s) { mApexPos=s; return this; }
  public LcMsFeature setPointWidth(int s) { mPointWidth=s; return this; }
  public LcMsPoint setAbundance(double s) { mAbundance=s; return this; }
  public LcMsFeature setArea(double s) { mArea=s; return this; }
  public LcMsPoint setApex(double s) { mApex=s; return this; }
//  public LcMsPoint setMz(double s) { mMz=s; return this; }
  public LcMsPoint setMzStdev(double s) { mMzStdev=s; return this; }
  public LcMsPoint setMzStdevEx(double s) { mMzStdevEx=s; return this; }
  public LcMsPoint setExclusivity(double s) { mExclusivity=s; return this; }

  public LcMsFeature setInitialCentroidRt(double s) { mInitialCentroidRt=s; return this; }
  public LcMsFeature wasBasedOn(criteria s) { mPeakSelection=s; return this; }
  public LcMsFeature setApexRT(double s) { mApexRT=s; return this; }
  public LcMsFeature setSNR(double s) { mSNR=s; return this; }
  public LcMsFeature setSimilarity(double s) { mSimilarity=s; return this; }
  public LcMsFeature setLower(double s) { mLower=s; return this; }
  public LcMsFeature setUpper(double s) { mUpper=s; return this; }

  public LcMsFeature inMS1(boolean s) { mInMS1=s; return this; }

  public boolean isBetterThan(LcMsFeature selected, float apexR, float dpD, int presentD)
  {
    // has to be better all around
    return (getY()>selected.getY()*apexR && getSimilarity()>selected.getSimilarity()+dpD && isPresentIn()>selected.isPresentIn()+presentD);
  }
  public static Double sumArea(Collection<LcMsFeature> data)
  {
    if (data==null) return null;

    double sum = 0d;
    for (LcMsFeature xy : data) if (xy!=null) sum += xy.getArea();

    return sum;
  }
  @Override
  public String toString()
  {
    String out = "rt "+ Tools.d2s(getRT(), 2) + ", m/z " + Tools.d2s(getMz(), 4) + ", ai " +
        Tools.d2s(getIntensity(), 1) + ", msec " + Tools.d2s(getFillTime(), 1) + ", area "+
        Tools.d2s(getArea(),1) + ", snr "+Tools.d2s(getSNR(), 2)+", dp "+Tools.d2s(getSimilarity(),2)+", %ex "+
        Tools.d2s(getExclusivity(),2)+(inMS1()?", ms1":"");
    return out;
  }
}
