package org.ms2ms.alg;

import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.data.ms.ClusteringSettings;
import org.ms2ms.data.ms.LcMsMsDataset;
import org.ms2ms.io.MsIO;
import org.ms2ms.utils.Tools;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Collection;

/** Copyright 2014-2015 ms2ms.org
 *
 *  Algorithms for spectral clustering
 *  Author: wyu
 */
public class Clustering
{
  /** Locate a set of potential spectral clusters based solely on precursor and retention time. To be called after
   *  initial survey so the mz and RT index are in place
   *
   * @param data contains the definition of a set of LC-MS/MS runs from a given experiment
   * @return
   */
  public static LcMsMsDataset seekNuclei(LcMsMsDataset data, ClusteringSettings settings)
  {
    // check for the data requirement
    if (data==null || !Tools.isSet(data.getMzRtFileOffset()) || !Tools.isSet(data.getTicFileOffset())) return data;

    RandomAccessFile bin = null;
    try
    {
      // open the binary MSMS cache
      bin = data.getSpCacheFile(2);
      // starting from the most intense spectrum
      for (Double tic : data.getTicFileOffset().keySet())
      {
        // for each spectrum, cells a wide slice of spectra from the mz/rt index
        for (Long id : data.getTicFileOffset().get(tic))
        {
          MsnSpectrum      ms  = MsIO.read(bin, new MsnSpectrum());
          Collection<Long> ids = data.getFileoffsetsByMzRt(ms.getPrecursor().getMz(),Tools.front(ms.getRetentionTimes()).getTime());
        }
      }
      // close the cache when done
      if (bin!=null) bin.close();
    }
    catch (IOException ie)
    {

    }
    return data;
  }
}
