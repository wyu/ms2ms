package org.ms2ms.test.ms;

import com.google.common.collect.SortedSetMultimap;
import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.cluster.DenseSimilarityGraph;
import org.expasy.mzjava.core.ms.cluster.MSTClusterBuilder;
import org.expasy.mzjava.core.ms.cluster.SimEdge;
import org.expasy.mzjava.core.ms.cluster.SimilarityGraph;
import org.expasy.mzjava.core.ms.spectrasim.DpSimFunc;
import org.expasy.mzjava.core.ms.spectrasim.NdpSimFunc;
import org.expasy.mzjava.core.ms.spectrasim.SimFunc;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.junit.Assert;
import org.junit.Test;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.data.ms.MaxQuant;
import org.ms2ms.io.MsIO;
import org.ms2ms.r.Dataframe;
import org.ms2ms.test.TestAbstract;
import org.ms2ms.utils.Stats;
import org.ms2ms.utils.Tools;

import java.io.RandomAccessFile;
import java.util.Collection;
import java.util.List;
import java.util.Set;

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
  public void preclusterByMzRT() throws Exception
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
        MsIO.writeSpectra("/media/data/tmp/examples_"+spectra.size()+"_"+Tools.d2s(mz, 4)+"_"+Tools.d2s(rt, 2)+".ms2", spectra);
        System.out.println(spectra.size()+"...");
        if (++counts>5) return;
      }
    }
    System.out.println();
  }
//  @Test
//  public void Kmean() throws Exception
//  {
//    List<MsnSpectrum> spectra = MsIO.readSpectra("/tmp/examples_13_425.7649_41.43.ms2");
//
//    int minCountAbs=1; double minSim=0.5d, mzTol=0.05d, maxConsPkGrpMzWidth=0.05d, minCountPct=50d;
//    KMeansPeakListClusterer kmean = new KMeansPeakListClusterer(minSim, mzTol, maxConsPkGrpMzWidth, minCountAbs, minCountPct);
//
//    for (MsnSpectrum spec : spectra)
//    {
//      kmean.add(spec);
//    }
//    kmean.cluster();
//    System.out.println(kmean.getClusterMembers().size());
//  }
  public SimilarityGraph<MsnSpectrum> newSimGraph(Collection<MsnSpectrum> spectra)
  {
    if (!Tools.isSet(spectra)) return null;

    DenseSimilarityGraph.Builder<MsnSpectrum> builder = new DenseSimilarityGraph.Builder<MsnSpectrum>();
    SimFunc                                       sim = new NdpSimFunc(50, new AbsoluteTolerance(0.5));
    for (MsnSpectrum A : spectra)
    {
      for (MsnSpectrum B : spectra)
      {
        if (A!=B)
        {
          SimEdge<MsnSpectrum> edge = builder.add(A, B, sim.calcSimilarity(A, B));
          System.out.println(edge.getVertex1().getScanNumbers().toString()+" --> " + edge.getVertex2().getScanNumbers().toString() + ":" +Tools.d2s(edge.getScore(), 2));
        }
      }
    }

    return builder.build();
  }
  @Test
  public void mst() throws Exception
  {
    List<MsnSpectrum> spectra = MsIO.readSpectra("/media/data/tmp/examples_13_425.7649_41.43.ms2");
    SimilarityGraph<MsnSpectrum> graph = newSimGraph(spectra);
    MSTClusterBuilder clusterer = new MSTClusterBuilder(0.5);

    Collection<Set<String>> clusters = clusterer.cluster(graph);

    System.out.println(clusters.size());
    Assert.assertEquals(clusters.size(), 5);
//    int minCountAbs=1; double minSim=0.5d, mzTol=0.05d, maxConsPkGrpMzWidth=0.05d, minCountPct=50d;
  }
}
