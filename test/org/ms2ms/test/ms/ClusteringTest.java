package org.ms2ms.test.ms;

import com.google.common.collect.Range;
import com.google.common.collect.SortedSetMultimap;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.expasy.mzjava.proteomics.ms.cluster.KMeansPeakListClusterer;
import org.expasy.mzjava.proteomics.ms.cluster.MSTPeakListClusterer;
import org.junit.Test;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.data.ms.MaxQuant;
import org.ms2ms.io.MsIO;
import org.ms2ms.io.MsReaders;
import org.ms2ms.r.Dataframe;
import org.ms2ms.test.TestAbstract;
import org.ms2ms.utils.Stats;
import org.ms2ms.utils.Tools;

import java.io.RandomAccessFile;
import java.util.Collection;
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
  public void precluster() throws Exception
  {
    Dataframe msms = Dataframe.readtable("/media/data/maxquant/20081129_Orbi6_NaNa_SA_FASP_out/combined/txt/composite_scans.txt",'\t').setTitle("msms");

    MultiTreeTable<Double, Double, String> mz_rt_row = msms.index(MaxQuant.V_MZ,MaxQuant.V_RT);
    SortedSetMultimap<     Double, String> tic_row   = Tools.reverse(msms.index(MaxQuant.V_TIC));

    double mztol=0.05, rttol=1.5;
    RandomAccessFile bin = new RandomAccessFile("/media/data/test/mzXML/cache435790877685301.ms2", "r");
    // starting from the most intense scan
    int counts=0;
    for (Double t : tic_row.keySet())
    {
      for (String r : tic_row.get(t))
      {
        double mz = Stats.toDouble(msms.cell(r, MaxQuant.V_MZ)), rt = Stats.toDouble(msms.cell(r, MaxQuant.V_RT));
        Collection<String>  slice = mz_rt_row.subset(mz-mztol, mz+mztol, rt-rttol,rt+rttol);
        List<MsnSpectrum> spectra = MsIO.readSpectra(bin, msms.getLongCol(MaxQuant.V_OFFSET, slice));
        MsIO.writeSpectra("/tmp/examples_"+spectra.size()+"_"+Tools.d2s(mz, 4)+"_"+Tools.d2s(rt, 2)+".ms2", spectra);
        System.out.println(spectra.size()+"...");
        if (++counts>5) return;
      }
    }
    System.out.println();
  }

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
