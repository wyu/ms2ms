package org.ms2ms.test;

import org.apache.hadoop.hbase.HBaseConfiguration;
import org.apache.hadoop.hbase.client.*;
import org.apache.hadoop.hbase.util.Bytes;
import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.PpmTolerance;
import org.expasy.mzjava.core.ms.peaklist.DoublePeakList;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.Peak;
import org.expasy.mzjava.proteomics.ms.spectrum.LibrarySpectrum;
import org.expasy.mzjava.proteomics.ms.spectrum.PepLibPeakAnnotation;
import org.junit.Test;
import org.ms2ms.alg.Peaks;
import org.ms2ms.mimsl.MIMSL;
import org.ms2ms.mzjava.AnnotatedSpectrum;
import org.ms2ms.nosql.HBasePeakList;
import org.ms2ms.nosql.HBaseProteomics;
import org.ms2ms.splib.SpLibs;
import org.ms2ms.utils.Tools;

import java.io.*;
import java.util.*;

/** Reading the content of splib
 *
 * Created by wyu on 4/13/14.
 */
public class SpLibsTest extends TestAbstract
{
  @Test
  public void readMsp() throws IOException
  {
//    Collection<LibrarySpectrum> spectra = SpLibs.readMsp(new File("/media/data/splib/human_crp_consensus_final_true_lib.msp"));
//    Collection<LibrarySpectrum> spectra = SpLibs.readMsp(new File("/media/data/splib/nist_nci_stdmix_consensus_final_true_lib.msp"));
    Collection<LibrarySpectrum> spectra = SpLibs.readMsp(new File("/media/data/splib/NIST_human_IT_2012-05-30.msp"));

    // save the spectrum and indice to HBase
    HBaseProteomics.index(spectra, HBasePeakList.SPEC_TRAP_CID, 50d, 450d, 7, 4d);
    HBaseProteomics.listTables();

    //assert spectra.size()==92;
  }
  @Test
  public void prepareMsp() throws IOException
  {
    // ensure that the table has been created
    HBaseProteomics.ensureTables();

    String root = "/media/data/splib/2013";
    HBaseProteomics.prepareMsps(root, HBasePeakList.SPEC_TRAP_CID, 50d, 450d, 7, 4d, "bsa_consensus_final_true_lib.msp",
        "human_crp_consensus_final_true_lib.msp","c_elegans_consensus_final_true_lib.msp","chicken_consensus_final_true_lib.msp",
        "dradiodurans_consensus_final_true_lib.msp","msmegmatis_consensus_final_true_lib.msp","nist_nci_stdmix_consensus_final_true_lib.msp");
    HBaseProteomics.prepareMsps(root, HBasePeakList.SPEC_TRAP_CID, 50d, 450d, 7, 4d,
        "drosophila_consensus_final_true_lib.msp","e_coli_consensus_final_true_lib.msp","human_b2mg_consensus_final_true_lib.msp",
        "yeast_pombe_consensus_final_true_lib.msp","sigmaups1_consensus_final_true_lib.msp","hsa201244f.msp","human_consensus_final_true_lib.msp",
        "rat_consensus_final_true_lib.msp","yeast_consensus_final_true_lib.msp","mouse_consensus_final_true_lib.msp");
    HBaseProteomics.prepareMsps(root, HBasePeakList.SPEC_TRAP_CID, 50d, 450d, 7, 4d,
        "yeast_pombe_consensus_final_true_lib.msp","sigmaups1_consensus_final_true_lib.msp","hsa_consensus_final_true_lib.msp","human_consensus_final_true_lib.msp",
        "rat_consensus_final_true_lib.msp","yeast_consensus_final_true_lib.msp","mouse_consensus_final_true_lib.msp");
    HBaseProteomics.prepareMsps(root, HBasePeakList.SPEC_TRAP_CID, 50d, 450d, 7, 4d,
        "hsa_consensus_final_true_lib.msp","human_consensus_final_true_lib.msp",
        "rat_consensus_final_true_lib.msp","yeast_consensus_final_true_lib.msp","mouse_consensus_final_true_lib.msp");

    HBaseProteomics.prepareMsps(root, HBasePeakList.SPEC_QTOF, 50d, 450d, 7, 4d, "rat_qtof_consensus_final_true_lib.msp","yeast_qtof_consensus_final_true_lib.msp");

    // peak annotation format is invalid? Exception thrown
//    HBaseProteomics.prepareMsps(root, HBasePeakList.SPEC_TRAP_HCD, 50d, 450d, 7, 4d, "human_hcd_selected_final_true_lib.msp");
    HBaseProteomics.listTables();
  }
  @Test
  public void readAllMsMsIndex() throws IOException
  {
    HConnection conn = HConnectionManager.createConnection(HBaseConfiguration.create());
    // get the number of row. Can be very expansive for a large table!!
    HTableInterface table = conn.getTable(HBasePeakList.TBL_MSMSINDEX);

    Scan scan = new Scan();
    scan.setCaching(1);
    scan.setBatch(1);
    scan.addFamily(Bytes.toBytes(HBasePeakList.FAM_PROP));

    ResultScanner resultScanner = table.getScanner(scan);
    Iterator<Result> iterator = resultScanner.iterator();
    byte[] last_row=null;
    while (iterator.hasNext())
    {
      Result next = iterator.next();
      for(Map.Entry<byte[], NavigableMap<byte[], NavigableMap<Long, byte[]>>> columnFamilyMap : next.getMap().entrySet())
      {
        for (Map.Entry<byte[], NavigableMap<Long, byte[]>> entryVersion : columnFamilyMap.getValue().entrySet())
        {
          for (Map.Entry<Long, byte[]> entry : entryVersion.getValue().entrySet())
          {
            byte[] row = next.getRow();
            float prec_mz = Bytes.toFloat(Bytes.head(row, 4));
            int pre_z = (int )row[4];
            String column = Bytes.toString(entryVersion.getKey());
            byte[] value = entry.getValue();
            long timesstamp = entry.getKey();
            if (last_row==null || last_row!=row) System.out.print("    row: " + prec_mz + ",+" + pre_z);
            System.out.print(";" + column);
            try
            {
              System.out.print("=" + Tools.d2s(Bytes.toDouble(value), 3));
            }
            catch (IllegalArgumentException e)
            {
              System.out.print("=" + Bytes.toInt(value));
            }
            if (last_row==null || last_row!=row) System.out.println();
            last_row=row;
          }
        }
      }
    }
  }
  @Test
  public void readAllSpectra() throws IOException
  {
    Collection<AnnotatedSpectrum> spectra = HBaseProteomics.query(Peaks.newPeakList(500.5d, 2), HBasePeakList.SPEC_TRAP_CID, new AbsoluteTolerance(0.5), 0d);
    assert spectra.size()>0;
  }
  @Test
  public void queryTest() throws IOException
  {
    // 500.730, +2(6): 318.20,520.19,568.30,683.25,782.32,869.35,
    PeakList<PepLibPeakAnnotation> ions = Peaks.newPeakList(500.73d, 2, "318.2", "568.3"); //782.32,869.35
    List<AnnotatedSpectrum>  candidates = MIMSL.run(ions, HBasePeakList.SPEC_TRAP_CID, new PpmTolerance(500d), new AbsoluteTolerance(0.5));

    System.out.println(MIMSL.printCandidates(null, candidates));
  }
}
