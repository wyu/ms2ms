package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.algo.Peaks;
import org.ms2ms.utils.Strs;

public class Ms2Pointer implements Comparable<Ms2Pointer>
{
  public String     run;
  public int        scan, z, hcode, npks, npks_upper;
  public float      mz, rt, prob, dp, mz_off;
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
    mz  = (float )ms.getPrecursor().getMz(); mz_off=0f;
    z   = ms.getPrecursor().getCharge();
    rt  = (float )ms.getRetentionTimes().getFirst().getTime()/60f;

    hcode=run.hashCode()+scan+Float.hashCode(mz);
  }
  public Ms2Pointer setMzOffset(float s) { mz_off=s; return this; }
  public float getMH() { return Peaks.toMH(mz,z); }

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

  public Ms2Pointer clone()
  {
    Ms2Pointer cloned = new Ms2Pointer(run, scan);

    cloned.z=z; cloned.npks=npks; cloned.npks_upper=npks_upper;
    cloned.mz=mz; cloned.rt=rt; cloned.prob=prob; cloned.dp=dp; cloned.mz_off=mz_off;
    cloned.cluster=cluster;

    cloned.pointer=pointer; // to the reading position in a binary file

    return cloned;
  }
}
