package org.ms2ms.apps;

import org.apache.hadoop.hbase.client.HConnection;
import org.apache.hadoop.hbase.client.HConnectionManager;
import org.apache.hadoop.hbase.client.HTableInterface;
import org.ms2ms.nosql.HBases;
import org.ms2ms.nosql.ms.HBasePeakList;
import org.ms2ms.nosql.ms.HBaseProteomics;
import org.ms2ms.nosql.ms.HBaseSpLib;
import org.ms2ms.utils.TabFile;
import org.ms2ms.utils.Tools;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 6/11/14
 * Time: 11:09 PM
 * To change this template use File | Settings | File Templates.
 */
public class SpLibBuilder extends Apps
{
  String mLibSource, mLibFormat, mLibOrganism, mLibVersion, mLibType, mLibName, mCfg, mLibList, mLibRoot;

  public SpLibBuilder() { super(); mAppName = "SpLibBuilder"; mVersion = "v0.001"; }
  protected void processCommandLine(String args[])
  {
    System.out.println(mAppName);
    System.out.println("Build Date: Jul 3, 2014; 5:45pm");
    super.processCommandLine(args);

    for (int i = 0; i < args.length; i++)
    {
      if (args[i].length() == 0) continue;
//      HBaseSpLib.set(table, "NIST", "msp", HBases.HUMAN, "JUN2013", HBasePeakList.SPEC_TRAP_HCD, "human_hcd_selected_final_true_lib");
//      mLibSource   = option(mLibSource, "Source of the spectral library such as NIST", args, i, "-s", "-source");
//      mLibFormat   = option(mLibFormat, "Format such as MSP", args, i, "-f", "-format");
//      mLibOrganism = option(mLibOrganism, "Organism", args, i, "-g", "-organism");
//      mLibVersion  = option(mLibVersion, "Version or the date when the library was constructed", args, i, "-v", "-version");
//      mLibType     = option(mLibType, "Spectrum Type such as CID/HCD", args, i, "-t", "-type");
//      mLibName     = option(mLibName, "Name of the library", args, i, "-n", "-name");
      mLibList     = option(mLibList, "A comma-delimited file containing a list of spectral libraries", args, i, "-L", "-libs");
      mLibRoot     = option(mLibRoot, "Root folder where the libraries are stored", args, i, "-in");
    }
  }
  @Override
  protected boolean doRun() throws Exception
  {
    if (!Tools.isSet(mLibList)) return false;

    HBases.verifyConnection();
    HBaseProteomics.ensureTables();

    HConnection      conn = HConnectionManager.createConnection(HBases.conf);
    HTableInterface table = conn.getTable(HBaseSpLib.TBL_SPLIB);

    TabFile cfg = new TabFile(mLibList, TabFile.comma);
    while (cfg.hasNext())
    {
      byte[] stype = HBasePeakList.SPEC_TRAP_CID;
      if      (cfg.get("spectype").equals("HCD" )) stype = HBasePeakList.SPEC_TRAP_HCD;
      else if (cfg.get("spectype").equals("QTOF")) stype = HBasePeakList.SPEC_QTOF;
      HBaseSpLib lib = HBaseSpLib.set(table, cfg.get("source"), cfg.get("format"), cfg.get("organism"), cfg.get("version"), stype, cfg.get("name"));
      System.out.println("\n" + HBaseProteomics.prepareLib(mLibRoot, lib, 50d, 450d, 7, 4d) + " entries prepared");
    }
    table.close(); conn.close();

    return true;
  }
  public static void main(String args[]) throws Exception
  {
    System.exit((new SpLibBuilder()).tryRun(args));
    // If we exit on the mac, the window goes away.
  }
}
