package org.ms2ms.apps;

import org.ms2ms.data.ms.LcMsMsDataset;
import org.ms2ms.data.ms.MaxQuant;
import org.ms2ms.utils.Tools;

/**
 * ** Copyright 2014-2015 ms2ms.org
 * <p/>
 * The initial step in the spectral clustering pipeline. It will survey the raw data files and setup the binary dump
 * and indices for the subsequent analysis
 * <p/>
 * Author: wyu
 * Date:   7/15/15
 */
public class Prep extends App
{
  private LcMsMsDataset lcmsms;

  @Override
  protected boolean doRun() throws Exception
  {
    lcmsms.init();
//    Dataframe msms = lcmsms.readMsMsWithAnnotations(),
//        survey = Dataframe.readtable(lcmsms.getResultRoot()+"scan_survey.txt",    '\t').setTitle("survey"),
//        offsets = Dataframe.merge(msms, survey, true, true, "Raw file", "Scan number").setTitle("offsets");

//    IOs.write(lcmsms.getResultRoot()+"composite_scans.txt", offsets.display().toString());

    return true;
  }
  protected void processCommandLine(String args[])
  {
    System.out.println(mAppName);
    System.out.println("Build Date: " + mBuild);
    super.processCommandLine(args);

    // let's build the main object first
    lcmsms=null;
    for (int i = 0; i < args.length; i++)
    {
      if (args[i].length()>0 && (args[i].equals("-x") || args[i].equals("-exec")))
      {
        // the source of the data
        if (Tools.equalsCaseless(args[i + 1], "mq")) lcmsms = new MaxQuant();
      }
    }
    if (lcmsms==null) throw new RuntimeException("LC-MS/MS processor not defined");

    for (int i = 0; i < args.length; i++)
    {
      if (args[i].length() == 0 || args[i].charAt(0) != '-') continue;
      // Is the argument a command-line option (starts with a '-') ?
      if      (args[i].equals("-r") || args[i].equals("-result"))
      {
        lcmsms.setResultsRoot(args[i+1]);
      }
      else if (args[i].equals("-w") || args[i].equals("-raw"))
      {
        lcmsms.setRawFileRoot(args[i+1]);
      }
    }
  }
  public static void main(String args[]) throws Exception
  {
    System.exit((new Prep()).tryRun(args));
    // If we exit on the mac, the window goes away.
  }
}
