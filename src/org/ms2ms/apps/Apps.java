package org.ms2ms.apps;

import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 6/11/14
 * Time: 11:06 PM
 * To change this template use File | Settings | File Templates.
 */
abstract public class Apps
{
  protected String                   mOutfileRoot;
  protected boolean                  mVerbose;
  protected List<String>             mUsages;
  //**************************************************************************
  // PRIVATE FIELDS
  //***-*********-*********-*********-*********-*********-*********-*********-
  protected String mAppName = "unTitled", mVersion;

  //--------------------------------------------------------------------------
  public int tryRun(String args[]) throws Exception
  {
    int exitStatus = 0;

    try
    {
      run(args);
    }
    catch (Throwable e)
    {
      e.printStackTrace();
      exitStatus = 1;
    }

    return(exitStatus);
  }

  public void run(String[] inArgs) throws Exception
  {
    processCommandLine(inArgs);

    Long msec0 = System.currentTimeMillis();
    // the actual process
    doRun();

    System.out.println("Duration (sec):" + Tools.d2s((System.currentTimeMillis() - msec0) * 0.001d, 2));
  }

  private boolean verbose() { return mVerbose; }

  protected void processCommandLine(String args[])
  {
    for (int i = 0; i < args.length; i++)
    {
      if (args[i].length() == 0) continue;

      // Is the argument a command-line option (starts with a '-') ?
      if (args[i].charAt(0) == '-')
      {
        if      (args[i].equals("-h") || args[i].equals("-help"))
        {
          usage();
          System.exit(0);
        }
        else if (args[i].equals("-d") || args[i].equals("-debug"))
        {
          setVerbose(true); // Turn verbose mode on
        }
        else if (args[i].equals("-v") || args[i].equals("-verbose"))
        {
          setVerbose(true); // Turn verbose mode on
        }
        else if (args[i].equals("-O") || args[i].equals("-outfile"))
        {
          mOutfileRoot = args[i+1];
        }
      }
    }
  }//--------------------------------------------------------------------------
  private void verboseMsg(String inMsg)
  {
    if (verbose())
    {
      System.out.println(inMsg);
    }
  }

  public void setVerbose(boolean inValue)
  {
    mVerbose = inValue;
  }

  protected void usage()
  {
    System.out.println("");
    System.out.println("USAGE FOR " + mAppName + ": " + mAppName + " [options]");
    System.out.println("\n    COMMAND-LINE OPTIONS:");

    if (Tools.isSet(mUsages))
      for (String u : mUsages)
        System.out.println("    " + u);
  }
  abstract protected boolean doRun() throws Exception;

  protected void close() throws IOException
  {
  }

  public String getAppName()         { return mAppName; }
  public void   setAppName(String s) { mAppName = s; }
  protected String option(String var, String usage, String[] args, int i, String... tags)
  {
    if (args[i].charAt(0)=='-' && Tools.isSet(args) && args.length>i+1 && Tools.isA(args[i], tags))
    {
      var=args[i+1];
      if (mUsages==null) mUsages = new ArrayList<String>();
      mUsages.add(Strs.toString(tags, "\t") + (Tools.isSet(usage) ? "\t:\t" + usage : ""));
    }
    return var;
  }
}