package org.ms2ms.test.ms;

import com.google.common.collect.Multimap;
import com.google.common.collect.SortedSetMultimap;
import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.cluster.*;
import org.expasy.mzjava.core.ms.spectrasim.NdpSimFunc;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.junit.Assert;
import org.junit.Test;
import org.ms2ms.algo.Clustering;
import org.ms2ms.algo.MsStats;
import org.ms2ms.algo.Spectra;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.data.ms.ClusteringSettings;
import org.ms2ms.data.ms.MaxQuant;
import org.ms2ms.data.ms.MsSettings;
import org.ms2ms.io.MsIO;
import org.ms2ms.r.Dataframe;
import org.ms2ms.test.TestAbstract;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Tools;

import java.io.RandomAccessFile;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: hliu
 * Date: 8/22/14
 * Time: 3:24 AM
 * To change this template use File | Settings | File Templates.
 */
public class ClusteringTest extends TestAbstract
{
  String root = "/Users/hliu/Desktop/App/2014/data/mzXML-centroid/";

  @Test
  public void clusterQC() throws Exception
  {
    Dataframe msms = abbr("/media/data/test/mzXML/composite_scans_clusters.txt");

    // check the quality of the cluster using the key IDs
  }
  @Test
  public void clustering() throws Exception
  {
    Dataframe msms = Clustering.cluster(abbr("/media/data/maxquant/20081129_Orbi6_NaNa_SA_FASP_out/combined/txt/composite_scans.txt"),
        new RandomAccessFile("/media/data/test/mzXML/cache173190685179316.ms2", "r"),
        new GlobalThresholdClusterBuilder<MsnSpectrum>(0.5),
        new NdpSimFunc(50, new AbsoluteTolerance(0.5)), new ClusteringSettings(MsSettings.ORBITRAP).setRt(1.5f, 3f));

    IOs.write("/media/data/test/mzXML/composite_scans_clusters.txt", msms.display("\t", "").toString());

/*
    // we're going to attempt total clustering. Just the first pass without optimization
    String[] selected = { MaxQuant.V_MZ,MaxQuant.V_RT,MaxQuant.V_TIC,MaxQuant.V_MSEQ,MaxQuant.V_SEQ,MaxQuant.V_OFFSET,MaxQuant.V_RAWFILE,MaxQuant.V_SCAN };
    Dataframe msms = Dataframe.readtable("/media/data/maxquant/20081129_Orbi6_NaNa_SA_FASP_out/combined/txt/composite_scans.txt", selected,'\t', true).setTitle("msms");

    MultiTreeTable<Double, Double, String> mz_rt_row = msms.index(MaxQuant.V_MZ,MaxQuant.V_RT);
    SortedSetMultimap<     Double, String> tic_row   = Tools.reverse(msms.index(MaxQuant.V_TIC));

    double mztol=0.05, rttol=1.5;
    Set<String> plural = new HashSet<>();
    RandomAccessFile bin = new RandomAccessFile("/media/data/test/mzXML/cache173190685179316.ms2", "r");
    GlobalThresholdClusterBuilder<MsnSpectrum> clusterer = new GlobalThresholdClusterBuilder<MsnSpectrum>(0.5);

    System.out.println("\nClustering the total MS/MS collection: " + msms.size());

    // starting from the most intense scan
    int counts=0, cluster_id=0;
    for (Double t : tic_row.keySet())
    {
      if (++counts%100  ==0) System.out.print(".");
      if (  counts%10000==0) System.out.println(cluster_id+"@"+Tools.d2s(100d*counts/tic_row.keySet().size(), 1)+"%");
      for (String r : tic_row.get(t))
      {
        if (plural.contains(r)) continue;

        double mz = MsStats.toDouble(msms.cell(r, MaxQuant.V_MZ)), rt = MsStats.toDouble(msms.cell(r, MaxQuant.V_RT));
        Collection<String>  slice = mz_rt_row.subset(mz-mztol, mz+mztol, rt-rttol,rt+rttol);
        // remove the rows already in the clusters (plural)
        Iterator<String> itr = slice.iterator();
        while (itr.hasNext())
        {
          if (plural.contains(itr.next())) itr.remove();
        }
        // compute the pair-wise sim scores and construct a simGraph. The rowid is stored as "Comment" with the spectrum
        SimilarityGraph<MsnSpectrum> graph = newSimGraph(MsIO.readSpectra(bin, msms.getLongColRow(MaxQuant.V_OFFSET, slice)));
        // form the clusters
        Collection<Set<MsnSpectrum>> clusters = clusterer.cluster(graph);
        // tag the spectra in cluster
        for (Set<MsnSpectrum> cluster : clusters)
        {
          if (cluster.size()<=1) continue;
          // start the cluster ID at 1
          cluster_id++;
          for (MsnSpectrum member : cluster)
          {
            // put it in the set so we can exclude them from further consideration
            plural.add(member.getComment());
            msms.put(member.getComment(), MaxQuant.V_CLUSTER, cluster_id);
          }
        }
      }
    }
    msms.init(Var.VarType.CATEGORICAL, MaxQuant.V_CLUSTER);
    System.out.println("Number of clusters: " + cluster_id);

    IOs.write("/media/data/test/mzXML/composite_scans_clusters.txt", msms.display("\t", "").toString());
*/
  }

  @Test
  public void preclusterByMSeq() throws Exception
  {
    String[] selected = { MaxQuant.V_MZ,MaxQuant.V_RT,MaxQuant.V_TIC,MaxQuant.V_MSEQ,MaxQuant.V_SEQ,MaxQuant.V_OFFSET,MaxQuant.V_Z };
    Dataframe msms = Dataframe.readtable("/media/data/maxquant/20081129_Orbi6_NaNa_SA_FASP_out/combined/txt/composite_scans.txt", selected,'\t', true).setTitle("msms");

    Multimap<String, String> seq_row   = msms.factorize(MaxQuant.V_MSEQ);

    RandomAccessFile bin = new RandomAccessFile("/media/data/test/mzXML/cache173190685179316.ms2", "r");
    // starting from the most intense scan
    int counts=0;
    for (String seq : seq_row.keySet())
    {
      if (seq_row.get(seq).size()<36) continue;

      Collection<String>  slice = seq_row.get(seq);
      Multimap<Integer, MsnSpectrum> z_spec = Spectra.toChargePeakList(MsIO.readSpectra(bin, msms.getLongCol(MaxQuant.V_OFFSET, slice)));
      for (Integer z : z_spec.keySet())
      {
        String fname = "/media/data/tmp/examples_"+z_spec.get(z).size()+"_"+seq+"_z"+z+".ms2";
        MsIO.writeSpectra(fname, z_spec.get(z));
        System.out.println("Writing to " + fname);
      }
      if (++counts>5) return;
    }
    System.out.println();
  }
  @Test
  public void preclusterByMzRT() throws Exception
  {
    String[] selected = { MaxQuant.V_MZ,MaxQuant.V_RT,MaxQuant.V_TIC,MaxQuant.V_MSEQ,MaxQuant.V_SEQ,MaxQuant.V_OFFSET };
    Dataframe msms = Dataframe.readtable("/media/data/maxquant/20081129_Orbi6_NaNa_SA_FASP_out/combined/txt/composite_scans.txt", selected,'\t', true).setTitle("msms");

    MultiTreeTable<Double, Double, String> mz_rt_row = msms.index(MaxQuant.V_MZ,MaxQuant.V_RT);
    SortedSetMultimap<     Double, String> tic_row   = Tools.reverse(msms.index(MaxQuant.V_TIC));

    double mztol=0.05, rttol=1.5;
    RandomAccessFile bin = new RandomAccessFile("/media/data/test/mzXML/cache173190685179316.ms2", "r");
    // starting from the most intense scan
    int counts=0;
    for (Double t : tic_row.keySet())
    {
      for (String r : tic_row.get(t))
      {
        double mz = MsStats.toDouble(msms.cell(r, MaxQuant.V_MZ)), rt = MsStats.toDouble(msms.cell(r, MaxQuant.V_RT));
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
  @Test
  public void mst() throws Exception
  {
    List<MsnSpectrum> spectra = MsIO.readSpectra("/media/data/tmp/examples_58__AAVPSGASTGIYEALELR__z2.ms2");
    SimilarityGraph<MsnSpectrum> graph = Clustering.newSimGraph(spectra,new NdpSimFunc(50, new AbsoluteTolerance(0.5)));
    MSTClusterBuilder clusterer = new MSTClusterBuilder(0.5);

    Collection<Set<String>> clusters = clusterer.cluster(graph);

    System.out.println(clusters.size());
    Assert.assertEquals(clusters.size(), 5);
//    int minCountAbs=1; double minSim=0.5d, mzTol=0.05d, maxConsPkGrpMzWidth=0.05d, minCountPct=50d;
  }
  @Test
  public void mstAllZ() throws Exception
  {
    List<MsnSpectrum> spectra = MsIO.readSpectra("/media/data/tmp/examples_16__AAELIANSLATAGDGLIELR__z2.ms2");
    spectra.addAll(MsIO.readSpectra("/media/data/tmp/examples_25__AAELIANSLATAGDGLIELR__z3.ms2"));

    SimilarityGraph<MsnSpectrum> graph = Clustering.newSimGraph(spectra,new NdpSimFunc(50, new AbsoluteTolerance(0.5)));

    MSTClusterBuilder<MsnSpectrum> clusterer = new MSTClusterBuilder<MsnSpectrum>(0.5);

    Collection<Set<MsnSpectrum>> clusters = clusterer.cluster(graph);

    System.out.println(clusters.size());
    Assert.assertEquals(clusters.size(), 5);
//    int minCountAbs=1; double minSim=0.5d, mzTol=0.05d, maxConsPkGrpMzWidth=0.05d, minCountPct=50d;
  }
  @Test
  public void thresholdAllZ() throws Exception
  {
    List<MsnSpectrum> spectra = MsIO.readSpectra("/media/data/tmp/examples_16__AAELIANSLATAGDGLIELR__z2.ms2");
    spectra.addAll(MsIO.readSpectra("/media/data/tmp/examples_25__AAELIANSLATAGDGLIELR__z3.ms2"));

    SimilarityGraph<MsnSpectrum> graph = Clustering.newSimGraph(spectra,new NdpSimFunc(50, new AbsoluteTolerance(0.5)));
    GlobalThresholdClusterBuilder clusterer = new GlobalThresholdClusterBuilder(0.5);
//    MSTClusterBuilder clusterer = new MSTClusterBuilder(0.5);

    Collection<Set<String>> clusters = clusterer.cluster(graph);

    System.out.println(clusters.size());
    Assert.assertEquals(clusters.size(), 5);
//    int minCountAbs=1; double minSim=0.5d, mzTol=0.05d, maxConsPkGrpMzWidth=0.05d, minCountPct=50d;
  }
  private Dataframe abbr(String fname) throws Exception
  {
    String[] selected = { MaxQuant.V_MZ,MaxQuant.V_RT,MaxQuant.V_TIC,MaxQuant.V_MSEQ,MaxQuant.V_SEQ,MaxQuant.V_OFFSET,MaxQuant.V_RAWFILE,MaxQuant.V_SCAN,MaxQuant.V_CLUSTER };
    Dataframe msms = Dataframe.readtable(fname, selected,'\t', true).setTitle("msms");

    return msms;
  }
}
