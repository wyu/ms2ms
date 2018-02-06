package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.utils.Strs;

public class Ms2Pointer implements Comparable<Ms2Pointer>
{
  public String     run;
  public int        scan, z, hcode;
  public float      mz, rt, prob, dp;
  public Ms2Cluster cluster;

  public long       pointer; // to the reading position in a binary file

  public Ms2Pointer(String r, int s)
  {
    super();
    run=r; scan=s;
  }
  public Ms2Pointer(String r, MsnSpectrum ms)
  {
    run =r;
    scan=ms.getScanNumbers().getFirst().getValue();
    mz  = (float )ms.getPrecursor().getMz();
    z   = ms.getPrecursor().getCharge();
    rt  = (float )ms.getRetentionTimes().getFirst().getTime()/60f;

    hcode=run.hashCode()+scan+Float.hashCode(mz);
  }

  @Override
  public int compareTo(Ms2Pointer o)
  {
    int c = run.compareTo(o.run);
    if (c==0) c=Integer.compare(scan, o.scan);
    if (c==0) c=  Float.compare(mz,o.mz);
    if (c==0) c=  Float.compare(rt,o.rt);

    return c;
  }

  @Override
  public int hashCode()
  {
    return hcode;
  }
  @Override
  public boolean equals(Object s)
  {
    if (s==null) return false;

    Ms2Pointer o = (Ms2Pointer )s;
    if (Strs.equals(run, o.run) && scan==o.scan && mz==o.mz) return true;
    return false;
  }
}
