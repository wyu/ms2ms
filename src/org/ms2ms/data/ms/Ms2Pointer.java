package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.algo.Peaks;
import org.ms2ms.data.Binary;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

public class Ms2Pointer implements Comparable<Ms2Pointer>, Binary, Ion
{
  public String     run, name;
  public int        scan, z, hcode, npks, npks_upper;
  public float      mz, rt, prob, dp, mz_off;
  public Ms2Cluster cluster;

  public long       pointer; // to the reading position in a binary file

  public Ms2Pointer() { super(); }
  public Ms2Pointer(String r, int s)
  {
    super();
    run=r; scan=s;
  }
  public Ms2Pointer(String r, MsnSpectrum ms)
  {
    run =r;
    scan=ms.getScanNumbers().size()>0?ms.getScanNumbers().getFirst().getValue():-1;
    mz  = (float )ms.getPrecursor().getMz(); mz_off=0f;
    z   = ms.getPrecursor().getCharge();
    rt  = Tools.isSet(ms.getRetentionTimes())?(float )ms.getRetentionTimes().getFirst().getTime()/60f:0f;

    if (Strs.isSet(ms.getComment()))    name = ms.getComment();
    if (Strs.isSet(ms.getFragMethod())) name = Strs.extend(name, ms.getFragMethod(), "$");

    hcode=run.hashCode()+scan+Float.hashCode(mz);
  }
  public Ms2Pointer setMzOffset(float s) { mz_off=s; return this; }
  public float getMH()     { return Peaks.toMH(mz,z); }
  public float getMz()     { return mz; }
  public int   getCharge() { return z; }

  @Override
  public int compareTo(Ms2Pointer o)
  {
    int c = (run!=null && o.run!=null) ? run.compareTo(o.run):0;
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

  @Override
  public void write(DataOutput ds) throws IOException
  {
    IOs.write(ds, run);
    IOs.write(ds, name);
    IOs.write(ds, scan);
    IOs.write(ds, z);
    IOs.write(ds, hcode);
    IOs.write(ds, npks);
    IOs.write(ds, npks_upper);
    IOs.write(ds, mz);
    IOs.write(ds, rt);
    IOs.write(ds, prob);
    IOs.write(ds, dp);
    IOs.write(ds, mz_off);
    IOs.write(ds, pointer);
  }

  @Override
  public void read(DataInput ds) throws IOException
  {
    run       =IOs.read(ds, run);
    name      =IOs.read(ds, name);
    scan      =IOs.read(ds, 0);
    z         =IOs.read(ds, 0);
    hcode     =IOs.read(ds, 0);
    npks      =IOs.read(ds, 0);
    npks_upper=IOs.read(ds, 0);
    mz        =IOs.read(ds, 0f);
    rt        =IOs.read(ds, 0f);
    prob      =IOs.read(ds, 0f);
    dp        =IOs.read(ds, 0f);
    mz_off    =IOs.read(ds, 0f);
    pointer   =IOs.read(ds, 0l);

    // cluster not persisted;
  }
  @Override
  public String toString()
  {
    return (cluster!=null?"$$":"")+(Strs.isSet(name)?(name+"::"):"")+run+(scan>0?"#"+scan:"")+"|z"+z+"|m/z"+ Tools.d2s(mz, 4)+"|min"+Tools.d2s(rt, 2)+
        (dp>0?"|dp"+Tools.d2s(dp,2):"")+(npks_upper>0?"|npks"+npks_upper:"");
  }
}
