package org.ms2ms.alg;

import com.google.common.collect.Multimap;
import com.google.common.collect.SortedSetMultimap;
import org.expasy.mzjava.core.ms.cluster.*;
import org.expasy.mzjava.core.ms.spectrasim.SimFunc;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.data.ms.ClusteringSettings;
import org.ms2ms.data.ms.LcMsMsDataset;
import org.ms2ms.data.ms.MaxQuant;
import org.ms2ms.io.MsIO;
import org.ms2ms.r.Dataframe;
import org.ms2ms.r.Var;
import org.ms2ms.utils.Tools;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.*;

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
  @Deprecated
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
  public static void QC(Dataframe data) throws Exception
  {
    // check the readiness of the incoming data
    if (data==null || !data.hasVars(MaxQuant.V_CLUSTER,MaxQuant.V_MSEQ)) return;

    // check the quality of the cluster using the peptide IDs
    Multimap<String, String> seq_row = data.factorize(MaxQuant.V_MSEQ);

    //
  }

  public static <T extends MsnSpectrum> Dataframe cluster(Dataframe msms, RandomAccessFile bin, ClusterBuilder<MsnSpectrum> clusterer, SimFunc sim, ClusteringSettings settings) throws Exception
  {
    // we're going to attempt total clustering. Just the first pass without optimization
//    String[] selected = { MaxQuant.V_MZ,MaxQuant.V_RT,MaxQuant.V_TIC,MaxQuant.V_MSEQ,MaxQuant.V_SEQ,MaxQuant.V_OFFSET,MaxQuant.V_RAWFILE,MaxQuant.V_SCAN };
//    Dataframe msms = Dataframe.readtable("/media/data/maxquant/20081129_Orbi6_NaNa_SA_FASP_out/combined/txt/composite_scans.txt", selected,'\t', true).setTitle("msms");

    MultiTreeTable<Double, Double, String> mz_rt_row = msms.index(MaxQuant.V_MZ,MaxQuant.V_RT);
    SortedSetMultimap<     Double, String> tic_row   = Tools.reverse(msms.index(MaxQuant.V_TIC));

//    double mztol=0.05, rttol=1.5;
    Set<String> plural = new HashSet<>();
//    RandomAccessFile bin = new RandomAccessFile("/media/data/test/mzXML/cache173190685179316.ms2", "r");
//    GlobalThresholdClusterBuilder<MsnSpectrum> clusterer = new GlobalThresholdClusterBuilder<MsnSpectrum>(0.5);

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

        Collection<String>  slice = mz_rt_row.subset(
            settings.getInstrument().getPrecursorTol().getMin(mz), settings.getInstrument().getPrecursorTol().getMax(mz),
            settings.getRtTol().getMin(rt), settings.getRtTol().getMax(rt));

        // remove the rows already in the clusters (plural)
        Iterator<String> itr = slice.iterator();
        while (itr.hasNext())
        {
          if (plural.contains(itr.next())) itr.remove();
        }
        // compute the pair-wise sim scores and construct a simGraph. The rowid is stored as "Comment" with the spectrum
        SimilarityGraph<MsnSpectrum> graph = newSimGraph(MsIO.readSpectra(bin, msms.getLongColRow(MaxQuant.V_OFFSET, slice)), sim);
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

    return msms;
  }

  /** An inefficient way to create a similarity graph.
   *
   * @param spectra are the incoming MS/MS spectra
   * @param sim defines the similarity. e.g. new NdpSimFunc(50, new AbsoluteTolerance(0.5))
   * @return the similarity graph to be send to a ClusterBuilder
   */
  public static <T extends MsnSpectrum> SimilarityGraph<T> newSimGraph(Collection<T> spectra, SimFunc sim)
  {
    if (!Tools.isSet(spectra)) return null;

    DenseSimilarityGraph.Builder<T> builder = new DenseSimilarityGraph.Builder<>();
    for (T A : spectra)
      for (T B : spectra)
        if (A!=B) builder.add(A, B, sim.calcSimilarity(A, B));

    return builder.build();
  }
}
