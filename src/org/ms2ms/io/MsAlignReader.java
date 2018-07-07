package org.ms2ms.io;

import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.expasy.mzjava.core.ms.spectrum.RetentionTimeDiscrete;
import org.expasy.mzjava.core.ms.spectrum.TimeUnit;
import org.ms2ms.data.NameValue;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.data.ms.LcMsMsInfo;
import org.ms2ms.mzjava.IsotopePeakAnnotation;
import org.ms2ms.r.Dataframe;
import org.ms2ms.utils.Strs;

import java.io.*;
import java.util.SortedSet;

public class MsAlignReader extends mzReader
{
  private BufferedReader mFileReader = null;
  private SortedSet<Double> mPrecursors;

  public MsAlignReader() { super(); }
  public MsAlignReader(String s)
  {
    try { mFileReader = new BufferedReader(new InputStreamReader(new FileInputStream(s))); }
    catch (Exception e) { e.printStackTrace();}
  }

  public BufferedReader getReader() { return mFileReader; }

  public MultiTreeTable<Double, Double, String> indexPeptidesFrMS2(Dataframe out, String filename)
  {
    MultiTreeTable<Double, Double, String> rt_mz_row = MultiTreeTable.create();
    for (String row : out.rows())
      if (out.cell(row,"RunFile").equals(filename) && out.getInteger(row, "MsLevel")>1)
        rt_mz_row.put(out.getDouble(row, "RT"), out.getDouble(row, "mz"), row);

    return rt_mz_row;
  }
  public MultiTreeTable<Double, Double, String> indexPeptidesFrPSM(Dataframe out, String mz_col, String RT_col)
  {
    MultiTreeTable<Double, Double, String> rt_mz_row = MultiTreeTable.create();
    for (String row : out.rows())
      rt_mz_row.put(out.getDouble(row, RT_col), out.getDouble(row, mz_col), row);

    return rt_mz_row;
  }
  public MsnSpectrum update(MsnSpectrum ms, LcMsMsInfo scans)
  {
    if (ms.getScanNumbers()==null || ms.getScanNumbers().size()==0)
      System.out.println();
    Integer scan = ms.getScanNumbers().getFirst().getValue();
    Double    rt = scans.getRT(null, scan);
    ms.addRetentionTime(new RetentionTimeDiscrete(rt*60d, TimeUnit.SECOND));

    ms.getPrecursor().setIntensity(scans.getAi(null, scan, 0d));
    return ms;
  }
  public Dataframe readPeptideFeatures(Dataframe out, LcMsMsInfo scans, double ppm, double dRT) throws IOException
  {
    MultiTreeTable<Double, Double, String> rt_mz_row = indexPeptidesFrPSM(out, "MH","RT");
    MultiTreeTable<Double, Double, Peak>   rt_mz_ms1 = MultiTreeTable.create(), rt_mz_mz1 = MultiTreeTable.create();

    // looping through the scans
    int rows=0;
    while (getReader().ready())
    {
      MsnSpectrum ms = readSpectrum(getReader());

      if (ms!=null && ms.size()>0 && update(ms, scans).getMsLevel()==1)
        XIC(ms, ppm, dRT, rt_mz_row, rt_mz_ms1, rt_mz_mz1);

      if (++rows%100==0) System.out.print(".");
    }
    System.out.println("$");

    // dump the traces for inspection
    FileWriter w = new FileWriter("/tmp/XICs.txt");
    w.write("rt0\tmz\trt\tai\n");
    writeXIC(w, rt_mz_ms1); w.close();

    // calculate the LC peak centroids
    for (String row : out.rows())
        out = detectXIC(out, row, "mz", "RT", rt_mz_ms1, rt_mz_mz1);

    return out;
  }
  public static MsnSpectrum readSpectrum(BufferedReader reader) throws IOException
  {
    // look for the beginning of the block
    while (reader.ready())
      if ("BEGIN IONS".equals(reader.readLine())) break;

    MsnSpectrum ms = new MsnSpectrum();
    Double mz=null; Integer z=null;
    while (reader.ready())
    {
      // look for the beginning of the block
      String line = reader.readLine();
      if ("END IONS".equals(line)) break; // quit at the end of the block
      if (line.indexOf("=")>0)
      {
        NameValue nv = NameValue.create(line, "=");
        if      ("SCANS".equals(nv.name)) ms.addScanNumber(nv.getInt());
        else if ("ACTIVATION".equals(nv.name)) ms.setFragMethod(nv.val);
        else if ("PRECURSOR_MZ".equals(nv.name)) mz=nv.getNumber();
        else if ("PRECURSOR_CHARGE".equals(nv.name)) z=nv.getInt();
      }
      else
      {
        String[] items = Strs.split(line, '\t');
        if (items.length>2)
          ms.add(Double.parseDouble(items[0]), Double.parseDouble(items[1]),
              new IsotopePeakAnnotation(Integer.parseInt(items[2]),0,0d));
      }
    }
    if (mz!=null && z!=null)
    {
      ms.setPrecursor(new Peak(mz,0d, z));
      ms.setMsLevel(2);
    }
    else ms.setMsLevel(1);

    return ms;

//    BEGIN IONS
//    ID=0
//    SCANS=1
//    444.1112529570313       926380.5027319421       1
//    518.1297466093749       464311.5385638939       1

//    BEGIN IONS
//    ID=711
//    SCANS=712
//    ACTIVATION=ETD
//    PRECURSOR_MZ=539.0160240992188
//    PRECURSOR_CHARGE=5
//    PRECURSOR_MASS=2690.0437404960944
//    1662.9145848671874      31613.025092928903      2
  }
}
