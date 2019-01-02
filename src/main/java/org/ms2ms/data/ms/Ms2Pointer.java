package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.algo.Peaks;
import org.ms2ms.data.Binary;
import org.ms2ms.r.Dataframe;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.FileWriter;
import java.io.IOException;

public class Ms2Pointer implements Comparable<Ms2Pointer>, Binary, Ion
{
  public String     run, name, id;
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
    if (ms.size()>0)
      for (int i=0; i<ms.size(); i++)
        hcode+=Double.hashCode(ms.getMz(i)) + Double.hashCode(ms.getIntensity(i))+i;
  }
  public Ms2Pointer setMzOffset(float s) { mz_off=s; return this; }
  public float getMH()     { return Peaks.toMH(mz,z); }
  public float getMz()     { return mz; }
  public int   getCharge() { return z; }

  public Dataframe details(Dataframe df, String rowid)
  {
    if (Strs.isSet(run))  df.put(rowid,"Run", run);
    if (Strs.isSet(name)) df.put(rowid,"Name", name);
    if (Strs.isSet(name)) df.put(rowid,"Cell", Strs.split(name,'$')[0]);
    if (Strs.isSet(id))   df.put(rowid,"ID", id);
    if (scan          >0) df.put(rowid,"Scan", scan+"");
    if (z            !=0) df.put(rowid,"z",z+"");
    if (npks          >0) df.put(rowid,"npks",npks+"");
    if (mz            >0) df.put(rowid,"mz", Tools.d2s(mz,4));
    if (rt           !=0) df.put(rowid,"RT",Tools.d2s(rt,2));
    if (prob          >0) df.put(rowid,"Prob",Tools.d2s(prob,2));
    if (dp            >0) df.put(rowid,"dp",Tools.d2s(dp,2));
    if (mz_off       !=0) df.put(rowid,"mz_offset",Tools.d2s(mz_off,4));
    if (pointer      !=0) df.put(rowid,"Pointer",pointer+"");

    return df;
  }
  public void details(FileWriter w) throws IOException
  {
    w.write((Strs.isSet(run)?run:"")+"\t");
    w.write((Strs.isSet(name)?name:"")+"\t");
    w.write((Strs.isSet(name)?Strs.split(name,'$')[0]:"")+"\t");
    w.write((Strs.isSet(id)?id:"")+"\t");
    w.write(scan+"\t");
    w.write(z+"\t");
    w.write(npks+"\t");
    w.write(Tools.d2s(mz,4)+"\t");
    w.write(Tools.d2s(rt,2)+"\t");
    w.write(Tools.d2s(prob,2)+"\t");
    w.write(Tools.d2s(dp,2)+"\t");
    w.write(Tools.d2s(mz_off,4)+"\t");
    w.write(pointer+"\n");
  }

  @Override
  public int compareTo(Ms2Pointer o)
  {
    int c = (run!=null && o.run!=null) ? run.compareTo(o.run):0;

    // only check the following if they are valid
    if (scan>=0 && mz>0 && rt>0)
    {
      if (c==0) c=Integer.compare(scan, o.scan);
      if (c==0) c=  Float.compare(mz,o.mz);
      if (c==0) c=  Float.compare(rt,o.rt);
    }
    else
    {
      // arbitary order
      if (c==0) c=Integer.compare(npks,  o.npks);
      if (c==0) c=Integer.compare(hcode, o.hcode);
    }

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
    if (Strs.equals(run, o.run) && scan>=0 && mz>0 && scan==o.scan && mz==o.mz) return true;
    return hcode==o.hcode;
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
    return (cluster!=null?"$$":"")+(Strs.isSet(name)?(name+"::"):"")+run+(scan>0?"#"+scan:"")+
        (z!=0?"|z"+z:"")+(mz!=0?"|m/z"+ Tools.d2s(mz, 4):"")+
        (rt>0?"|min"+Tools.d2s(rt, 2):"")+
        (dp>0?"|dp"+Tools.d2s(dp,2):"")+(npks_upper>0?"|npks"+npks_upper:"");
  }
}
