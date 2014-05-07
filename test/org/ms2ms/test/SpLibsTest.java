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
import org.ms2ms.mimsl.MIMSL;
import org.ms2ms.mzjava.AnnotatedSpectrum;
import org.ms2ms.nosql.HBasePeakList;
import org.ms2ms.nosql.HBaseProteomics;
import org.ms2ms.splib.SpLibs;
import org.ms2ms.utils.Tools;

import java.io.*;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.NavigableMap;

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
    HBaseProteomics.index(spectra, 50d, 450d, 7, 4d);
    HBaseProteomics.listTables();

    //assert spectra.size()==92;
  }
  @Test
  public void prepareMsp() throws IOException
  {
//    HBaseProteomics.prepareMsps("/media/data/splib", 50d, 450d, 7, 4d, "human_crp_consensus_final_true_lib.msp", "nist_nci_stdmix_consensus_final_true_lib.msp", "NIST_human_IT_2012-05-30.msp");
    HBaseProteomics.prepareMsps("/media/data/splib", 50d, 450d, 7, 4d, "NIST_mouse_IT_2012-04-21.msp","NIST_rat_IT_2012-04-16.msp","NIST_sigmaups1_IT_2011-05-24.msp","NIST_yeast_IT_2012-04-06.msp");
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
    Collection<AnnotatedSpectrum> spectra = HBaseProteomics.query(newPeakList(500.5d, 2), new AbsoluteTolerance(0.5), 0d);
    assert spectra.size()>0;
  }
  @Test
  public void queryTest() throws IOException
  {
    // 500.730, +2(6): 318.20,520.19,568.30,683.25,782.32,869.35,
    PeakList<PepLibPeakAnnotation> ions = newPeakList(500.73d, 2, "318.2", "568.3");
    MIMSL.run(ions, new PpmTolerance(500d), new AbsoluteTolerance(0.5));
  }

  /** makeup a peaklist using string shorthand
   *
   * @param mz and z are the m/z value and charge state of the precursor.
   * @param frags variable number of fragment ions. e.g. "334.5", "562/23", "mz/ai". Only the mz is required
   * @return
   */
  private static PeakList<PepLibPeakAnnotation> newPeakList(double mz, int z, String... frags)
  {
    PeakList<PepLibPeakAnnotation> spec = new DoublePeakList<PepLibPeakAnnotation>();
    spec.setPrecursor(new Peak(mz, 0d, z));
    // go thro the fragment ions if any
    for (String f : frags)
    {
      try
      {
        String[] strs = f.split("/");
        spec.add(Double.valueOf(strs[0]), strs.length>2?Double.valueOf(strs[2]):0d);
      }
      catch (Exception e) { e.printStackTrace(); }
    }
    return spec;
  }
}
