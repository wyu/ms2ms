package org.ms2ms.test.ms;

import org.apache.hadoop.hbase.HBaseConfiguration;
import org.apache.hadoop.hbase.client.*;
import org.apache.hadoop.hbase.util.Bytes;
import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.proteomics.ms.consensus.PeptideConsensusSpectrum;
import org.junit.Test;
import org.ms2ms.alg.Peaks;
import org.ms2ms.mzjava.AnnotatedSpectrum;
import org.ms2ms.nosql.HBase;
import org.ms2ms.nosql.HBasePeakList;
import org.ms2ms.nosql.HBaseProteomics;
import org.ms2ms.nosql.HBaseSpLib;
import org.ms2ms.splib.SpLibs;
import org.ms2ms.test.TestAbstract;
import org.ms2ms.utils.Tools;

import java.io.*;
import java.net.URISyntaxException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/** Reading the content of splib
 *
 * Created by wyu on 4/13/14.
 */
public class SpLibsTest extends TestAbstract
{
  @Test
  public void readMsp() throws IOException
  {
//    Collection<PeptideConsensusSpectrum> spectra = SpLibs.readMsp(new File("/media/data/splib/human_crp_consensus_final_true_lib.msp"));
//    Collection<PeptideConsensusSpectrum> spectra = SpLibs.readMsp(new File("/media/data/splib/nist_nci_stdmix_consensus_final_true_lib.msp"));
    Collection<PeptideConsensusSpectrum> spectra = SpLibs.readMsp(new File("/media/data/splib/NIST_human_IT_2012-05-30.msp"));

    // save the spectrum and indice to HBase
    HBaseProteomics.index(spectra, HBasePeakList.SPEC_TRAP_CID, 50d, 450d, 7, 4d);
    HBaseProteomics.listTables();

    //assert spectra.size()==92;
  }
//  @Test
//  public void readSpLib() throws IOException
//  {
//    long n = HBaseProteomics.prepareLib(new File("/media/data/splib/2013/human/HumanPlasma_2012-08_all.sptxt"), HBasePeakList.SPEC_TRAP_CID, 20, 450, 0, 0);
//    assert n==92;
//  }
/*
  @Test
  public void prepareLibs() throws IOException
  {
    // ensure that the table contains been created
    HBaseProteomics.ensureTables();

    String root = "/media/data/splib/2013";
    HBaseProteomics.prepareLibs(root, HBasePeakList.SPEC_TRAP_CID, 50d, 450d, 7, 4d, "bsa_consensus_final_true_lib.msp",
        "human_crp_consensus_final_true_lib.msp", "c_elegans_consensus_final_true_lib.msp", "chicken_consensus_final_true_lib.msp",
        "dradiodurans_consensus_final_true_lib.msp", "msmegmatis_consensus_final_true_lib.msp", "nist_nci_stdmix_consensus_final_true_lib.msp");
    HBaseProteomics.prepareLibs(root, HBasePeakList.SPEC_TRAP_CID, 50d, 450d, 7, 4d,
        "drosophila_consensus_final_true_lib.msp", "e_coli_consensus_final_true_lib.msp", "human_b2mg_consensus_final_true_lib.msp",
        "yeast_pombe_consensus_final_true_lib.msp", "sigmaups1_consensus_final_true_lib.msp", "hsa201244f.msp", "human_consensus_final_true_lib.msp",
        "rat_consensus_final_true_lib.msp", "yeast_consensus_final_true_lib.msp", "mouse_consensus_final_true_lib.msp");
    HBaseProteomics.prepareLibs(root, HBasePeakList.SPEC_TRAP_CID, 50d, 450d, 7, 4d,
        "yeast_pombe_consensus_final_true_lib.msp", "sigmaups1_consensus_final_true_lib.msp", "hsa_consensus_final_true_lib.msp", "human_consensus_final_true_lib.msp",
        "rat_consensus_final_true_lib.msp", "yeast_consensus_final_true_lib.msp", "mouse_consensus_final_true_lib.msp");
    HBaseProteomics.prepareLibs(root, HBasePeakList.SPEC_TRAP_CID, 50d, 450d, 7, 4d,
        "hsa_consensus_final_true_lib.msp", "human_consensus_final_true_lib.msp",
        "rat_consensus_final_true_lib.msp", "yeast_consensus_final_true_lib.msp", "mouse_consensus_final_true_lib.msp");

    HBaseProteomics.prepareLibs(root, HBasePeakList.SPEC_QTOF, 50d, 450d, 7, 4d, "rat_qtof_consensus_final_true_lib.msp", "yeast_qtof_consensus_final_true_lib.msp");

    // peak annotation format is invalid? Exception thrown
//    HBaseProteomics.prepareMsps(root, HBasePeakList.SPEC_TRAP_HCD, 50d, 450d, 7, 4d, "human_hcd_selected_final_true_lib.msp");
    HBaseProteomics.listTables();
  }
*/
  @Test
  public void prepareAllLibs() throws IOException, URISyntaxException
  {
    HBaseProteomics.ensureTables();

    Collection<HBaseSpLib> libs = HBaseSpLib.getAll();
    for (HBaseSpLib lib : libs)
    {
      // for debugging only
      if (lib==null || lib.getName().indexOf("HumanPlasma_2012-08_all")==0) continue;
      if (lib==null || lib.getName().indexOf("nist_nci_stdmix_consensus_final_true_lib")<0) continue;
      // preparing the library
      System.out.println("\n" + HBaseProteomics.prepareLib("/media/data/splib/2013", lib, 50d, 450d, 7, 4d) + " entries prepared");
    }
/*
    HBaseProteomics.ensureTables();

    HConnection      conn = HConnectionManager.createConnection(HBaseConfiguration.create());
    HTableInterface table = conn.getTable(HBaseSpLib.TBL_SPLIB);
    Scan             scan = new Scan();

    scan.setCaching(1); scan.setBatch(1);
    scan.addColumn(HBase.FAM_ID, HBase.COL_ENDID);

    int rows = 0;
    ResultScanner resultScanner = table.getScanner(scan);
    Iterator<Result>   iterator = resultScanner.iterator();
    while (iterator.hasNext())
    {
      Result result = iterator.next();
      //if you want to cells the entire row
      Get cells = new Get(result.row());
      Result entireRow = table.cells(cells);
      // process the row
      HBaseSpLib lib = new HBaseSpLib(entireRow);
      System.out.println("row " + (rows++) + ": " + lib.toString());
      // for debugging only
      if (lib==null || lib.getName().indexOf("HumanPlasma_2012-08_all")<0) continue;
      // preparing the library
      long counts = HBaseProteomics.prepareLib("/media/data/splib/2013", lib, 50d, 450d, 7, 4d);
    }
    table.close(); conn.close();
*/
  }
  @Test
  public void readAllMsMsIndex() throws IOException
  {
    HConnection conn = HConnectionManager.createConnection(HBaseConfiguration.create());
    // cells the number of row. Can be very expansive for a large table!!
    HTableInterface table = conn.getTable(HBase.TBL_MSMSINDEX);

    Scan scan = new Scan();
    scan.setCaching(1);
    scan.setBatch(1);
    scan.addFamily(HBase.FAM_PROP);

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
  public  void setSpLibs() throws IOException
  {
    HBaseProteomics.ensureTables();

    HConnection      conn = HConnectionManager.createConnection(HBaseConfiguration.create());
    HTableInterface table = conn.getTable(HBaseSpLib.TBL_SPLIB);

    HBaseSpLib.set(table, "PeptideAtlas", "sptxt", HBase.HUMAN, "AUG2012", HBasePeakList.SPEC_TRAP_CID, "HumanPlasma_2012-08_all");
    HBaseSpLib.set(table, "PeptideAtlas", "sptxt", HBase.HUMAN, "JUL2013", HBasePeakList.SPEC_TRAP_CID, "ISB_Hs_phospho_20130715");
    HBaseSpLib.set(table, "PeptideAtlas", "sptxt", HBase.HUMAN, "JUL2013", HBasePeakList.SPEC_TRAP_CID, "ISB_Hs_phospho_SEMI_20130715");
    HBaseSpLib.set(table, "PeptideAtlas", "sptxt", HBase.YEAST, "JUL2013", HBasePeakList.SPEC_TRAP_CID, "ISB_Sc_Phospho_20130715");
    HBaseSpLib.set(table, "PeptideAtlas", "sptxt", HBase.MM, "JUL2013", HBasePeakList.SPEC_TRAP_CID, "ISB_Mm_Phospho_20130715");
    HBaseSpLib.set(table, "PeptideAtlas", "sptxt", HBase.WORM, "JUL2013", HBasePeakList.SPEC_TRAP_CID, "ISB_Ce_Phospho_20130715_UniProt201306");
    HBaseSpLib.set(table, "PeptideAtlas", "sptxt", HBase.DM, "JUL2013", HBasePeakList.SPEC_TRAP_CID, "ISB_Dm_Phospho_20130715_UniProt201306");

    HBaseSpLib.set(table, "NIST", "msp", HBase.HUMAN, "JUN2013", HBasePeakList.SPEC_TRAP_HCD, "human_hcd_selected_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.HUMAN, "JUN2010", HBasePeakList.SPEC_TRAP_CID, "hsa_consensus_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.HUMAN, "MAY2009", HBasePeakList.SPEC_TRAP_CID, "human_b2mg_consensus_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.HUMAN, "JUN2013", HBasePeakList.SPEC_TRAP_CID, "human_consensus_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.HUMAN, "MAY2011", HBasePeakList.SPEC_TRAP_CID, "human_crp_consensus_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.MOUSE, "MAY2013", HBasePeakList.SPEC_TRAP_CID, "mouse_consensus_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.RAT, "MAY2013", HBasePeakList.SPEC_TRAP_CID, "rat_consensus_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.RAT, "JUN2013", HBasePeakList.SPEC_QTOF, "rat_qtof_consensus_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.YEAST, "MAY2012", HBasePeakList.SPEC_TRAP_CID, "yeast_consensus_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.YEAST, "MAY2012", HBasePeakList.SPEC_QTOF, "yeast_qtof_consensus_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.POMBE, "MAY2013", HBasePeakList.SPEC_TRAP_CID, "yeast_pombe_consensus_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.ECOLI, "MAY2013", HBasePeakList.SPEC_TRAP_CID, "e_coli_consensus_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.BOVINE, "JUN2011", HBasePeakList.SPEC_TRAP_CID, "bsa_consensus_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.WORM, "MAY2011", HBasePeakList.SPEC_TRAP_CID, "c_elegans_consensus_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.CHICK, "MAY2011", HBasePeakList.SPEC_TRAP_CID, "chicken_consensus_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.MM, "JUN2010", HBasePeakList.SPEC_TRAP_CID, "msmegmatis_consensus_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.DRADI, "JUN2010", HBasePeakList.SPEC_TRAP_CID, "dradiodurans_consensus_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.DROSO, "JUN2013", HBasePeakList.SPEC_TRAP_CID, "drosophila_consensus_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.MIXED, "JUN2010", HBasePeakList.SPEC_TRAP_CID, "nist_nci_stdmix_consensus_final_true_lib");
    HBaseSpLib.set(table, "NIST", "msp", HBase.MIXED, "MAY2011", HBasePeakList.SPEC_TRAP_CID, "sigmaups1_consensus_final_true_lib");

    table.close(); conn.close();
  }
  @Test
  public void peakMatcher()
  {
    String deci = "110.0717	4108.5	\"?\"", cid = "373.3\t144\t\"?i 21/25 0.4\"", exp="115.0868\t1.91984e+006\t\"a2/1.7ppm,a4-28^2/1.7ppm,Int/AA-28/1.7ppm\"", line=exp;

    //Pattern peakPattern = Pattern.compile("^([+-]?\\d+\\.?\\d*(?:[eE][-+]?\\d+)?)\\s+([+-]?\\d+)\\s+\"([^\"]+)\"$");
    //Pattern peakPattern = Pattern.compile("^([+-]?\\d+\\.?\\d*(?:[eE][-+]?\\d+)?)\\s+([+-]?\\d+(?:\\.\\d+)?)\\s+\"([^\"]+)\"$");
    Pattern peakPattern   = Pattern.compile("^([+-]?\\d+\\.?\\d*(?:[eE][-+]?\\d+)?)\\s+([+-]?\\d+\\.?\\d*(?:[eE][-+]?\\d+)?)\\s+\"([^\"]+)\"$");

    Matcher matcher = peakPattern.matcher(line);

    if (matcher.matches())
    {
      System.out.println(matcher.group(1) + ", " + matcher.group(2) + ", " + matcher.group(3));
//      builder.addPeak(Double.parseDouble(matcher.group(1)), Double.parseDouble(matcher.group(2)), matcher.group(3));
    } else {

      throw new IllegalArgumentException(line + " is not a valid peak line");
    }
  }
}
