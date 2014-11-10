package org.ms2ms.test;

import com.google.common.collect.Range;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.expasy.mzjava.proteomics.ms.cluster.KMeansPeakListClusterer;
import org.expasy.mzjava.proteomics.ms.cluster.MSTPeakListClusterer;
import org.junit.Test;
import org.ms2ms.io.MsReaders;

import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: hliu
 * Date: 8/22/14
 * Time: 3:24 AM
 * To change this template use File | Settings | File Templates.
 */
public class ClusteringTest extends TestAbstract
{
  String root = "/Users/hliu/Desktop/Apps/2014/data/mzXML-centroid/";

  @Test
  public void prepare() throws Exception
  {
    List<MsnSpectrum> spectra = MsReaders.readMzXML(root + "20081129_Orbi6_NaNa_SA_FASP_blacktips_01.mzXML", Range.closed(20d, 21d), 2);

    int minCountAbs=1; double minSim=0.5d, mzTol=0.05d, maxConsPkGrpMzWidth=0.05d, minCountPct=50d;
    KMeansPeakListClusterer kmean = new KMeansPeakListClusterer(minSim, mzTol, maxConsPkGrpMzWidth, minCountAbs, minCountPct);
    MSTPeakListClusterer      mst = new MSTPeakListClusterer(   minSim, mzTol, maxConsPkGrpMzWidth, minCountAbs, minCountPct);

    for (MsnSpectrum spec : spectra)
    {
      kmean.add(spec);
    }
    kmean.cluster();
    System.out.println();
  }

  @Test
  public void Kmean() throws Exception
  {
    List<MsnSpectrum> spectra = MsReaders.readMzXML(root + "20081129_Orbi6_NaNa_SA_FASP_blacktips_01.mzXML", Range.closed(20d, 21d), 2);

    int minCountAbs=1; double minSim=0.5d, mzTol=0.05d, maxConsPkGrpMzWidth=0.05d, minCountPct=50d;
    KMeansPeakListClusterer kmean = new KMeansPeakListClusterer(minSim, mzTol, maxConsPkGrpMzWidth, minCountAbs, minCountPct);
    MSTPeakListClusterer      mst = new MSTPeakListClusterer(   minSim, mzTol, maxConsPkGrpMzWidth, minCountAbs, minCountPct);

    for (MsnSpectrum spec : spectra)
    {
      kmean.add(spec);
    }
    kmean.cluster();
    System.out.println();
  }
}
