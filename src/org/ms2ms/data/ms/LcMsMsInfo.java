package org.ms2ms.data.ms;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import org.ms2ms.r.Dataframe;
import org.ms2ms.utils.Tools;

public class LcMsMsInfo
{
  private Integer     mMinMsLevel=2;
  private Dataframe            scans;
  private Table<String, Integer, String> run_scan_row;

  public LcMsMsInfo() { super(); }
  public LcMsMsInfo(String s, Integer ms_level)
  {
    scans = Dataframe.create(s, '\t'); mMinMsLevel=ms_level;
    prepare();
  }

  private boolean hasRow(String run, Integer scan) { return run_scan_row!=null && run_scan_row.contains(run,scan); }
  private Object     get(String run, Integer scan, String col, Object _def)
  {
    if (run==null)
    {
      return (scan!=null && Tools.isSet(run_scan_row.column(scan))) ?
          scans.cell(Tools.front(run_scan_row.column(scan).values()), col) : _def;
    }
    return hasRow(run,scan) && scans.cell(run_scan_row.get(run,scan), col)!=null ?
        scans.cell(run_scan_row.get(run,scan), col):_def;
  }
  private Double getDouble(String run, Integer scan, String col)
  {
    try { return (Double )get(run,scan,col, null); } catch (Exception e) { return null; }
  }
  private Double getDouble(String run, Integer scan, String col, Double _def)
  {
    Double val = getDouble(run,scan,col);
    return val!=null?val:_def;
  }
  public Double getMz(String run, Integer scan, Double _def) { return getDouble(run,scan,"XIC.mz",_def); }
  public Double getRT(String run, Integer scan, Double _def) { return getDouble(run,scan,"RT",_def); }
  // Intensity is the nominal value from the MS2 header
  public Double getAi(String run, Integer scan, Double _def) { return getDouble(run,scan,"XIC.area",_def); }

  public Double getMz(String run, Integer scan) { return getDouble(run,scan,"XIC.mz"); }
  public Double getRT(String run, Integer scan) { return getDouble(run,scan,"RT"); }
  // Intensity is the nominal value from the MS2 header
  public Double getAi(String run, Integer scan) { return getDouble(run,scan,"XIC.area"); }

  public int    size() { return scans!=null?scans.size():0; }

  private LcMsMsInfo prepare()
  {
    run_scan_row = HashBasedTable.create();

    // prepare the scan-based index
    for (String row : scans.rows())
      if (((Long )scans.row(row).get("MsLevel"))>=mMinMsLevel)
        run_scan_row.put(scans.row(row).get("Run").toString(), ((Long )scans.row(row).get("Scan")).intValue(), row);

    return this;
  }
}
