package org.ms2ms.io;

import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.*;
import org.ms2ms.utils.Tools;

import java.io.*;

/** Responsible for writing out the MS objects in various formats
 * User: wyu
 * Date: 9/29/14
 */
public class MsWriters
{
  // BufferedWriter
  public static void write(DataOutput w, PeakList ms) throws IOException
  {
    w.writeDouble(ms.getPrecursor().getMz());
    w.writeDouble(ms.getPrecursor().getIntensity());
    w.writeInt(   ms.getPrecursor().getCharge());

    w.writeInt(   ms.size());
    if (ms.size()>0)
    {
      for (int i=0; i<ms.size(); i++)
      {
        w.writeDouble(ms.getMz(i));
        w.writeDouble(ms.getIntensity(i));
      }
    }
  }
  public static void write(DataOutput w, MsnSpectrum ms) throws IOException
  {
    write(w, (PeakList)ms);

    int counts = ms.getRetentionTimes()!=null&&ms.getRetentionTimes().size()>0?ms.getRetentionTimes().size():0;
    w.writeInt(counts);
    if (counts>0)
      for (RetentionTime rt : ms.getRetentionTimes()) w.writeDouble(rt.getTime());
  }
  public static PeakList read(DataInput w, PeakList ms) throws IOException
  {
    if (ms!=null)
    {
      ms.setPrecursor(new Peak(w.readDouble(), w.readDouble(), w.readInt()));

      int counts = w.readInt();
      if (counts>0)
        for (int i=0; i<counts; i++)
          ms.add(w.readDouble(), w.readDouble());
    }
    return ms;
  }
  public static MsnSpectrum read(DataInput w, MsnSpectrum ms) throws IOException
  {
    if (ms!=null)
    {
      ms = (MsnSpectrum )read(w, (PeakList)ms);

      int counts = w.readInt();
      if (counts>0)
        for (int i=0; i<counts; i++)
          ms.addRetentionTime(new RetentionTimeDiscrete(w.readDouble(), TimeUnit.SECOND));
    }
    return ms;
  }
}
