package org.ms2ms.apps;

import org.ms2ms.utils.Strs;

/**
 * ** Copyright 2014-2015 ms2ms.org
 * <p/>
 * A standalone application to calculate the matching probability curve
 * from a collection of annotated MS/MS spectra
 * <p/>
 * Author: wyu
 * Date:   6/9/15
 */
public class Guru extends App
{
  private String mIdRoot;
  public Guru() { super(); mAppName = "Guru"; mVersion = "v0.001"; mBuild="Jun 9, 2015"; }
  protected void processCommandLine(String args[])
  {
    System.out.println(mAppName);
    System.out.println("Build Date: " + mBuild);
    super.processCommandLine(args);

    for (int i = 0; i < args.length; i++)
    {
      if (args[i].length() == 0) continue;
      mIdRoot     = option(mIdRoot, "The folder containing the mzXML files and msms.txt from MaxQuant", args, i, "-I", "-ids");
    }
  }
  public String getOutFile() { return null; }

  @Override
  protected boolean doRun() throws Exception
  {
    if (!Strs.isSet(mIdRoot)) return false;


    return true;
  }

  @Override
  protected void addProperty(String... vals)
  {

  }

  public static void main(String args[]) throws Exception
  {
    System.exit((new Guru()).tryRun(args));
    // If we exit on the mac, the window goes away.
  }
}
