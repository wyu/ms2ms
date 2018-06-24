package org.ms2ms.data.ms;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import org.ms2ms.r.Dataframe;

public class LcMsMsInfo
{
  private Dataframe            scans;
  private Table<String, Integer, String> run_scan_row;

  public LcMsMsInfo() { super(); }
  public LcMsMsInfo(String s)
  {
    scans = Dataframe.create(s, '\t');
    prepare();
  }

  private boolean hasRow(String run, Integer scan) { return run_scan_row!=null && run_scan_row.contains(run,scan); }
  private Object     get(String run, Integer scan, String col, Object _def)
  {
    return hasRow(run,scan) && scans.cell(run_scan_row.get(run,scan), col)!=null ?
        scans.cell(run_scan_row.get(run,scan), col):_def;
  }
  public Double getRT(String run, Integer scan, Double _def) { return (Double )get(run,scan,"RT",_def); }
  // Intensity is the nominal value from the MS2 header
  public Double getAi(String run, Integer scan, Double _def) { return (Double )get(run,scan,"XIC.area",_def); }
  public int    size() { return scans!=null?scans.size():0; }

  private LcMsMsInfo prepare()
  {
    run_scan_row = HashBasedTable.create();

    // prepare the scan-based index
    for (String row : scans.rows())
      if (((Long )scans.row(row).get("MsLevel"))==2)
        run_scan_row.put(scans.row(row).get("Run").toString(), ((Long )scans.row(row).get("Scan")).intValue(), row);

    return this;
  }
}
