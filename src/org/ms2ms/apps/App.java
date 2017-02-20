package org.ms2ms.apps;

import com.google.common.collect.BiMap;
import com.typesafe.config.Config;
import com.typesafe.config.ConfigFactory;
import org.ms2ms.graph.Property;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 6/11/14
 * Time: 11:06 PM
 * To change this template use File | Settings | File Templates.
 */
abstract public class App
{
  protected static BiMap<String, String> sParamKeys;

  protected Property                 mParameters;
  protected String                   mOutfileRoot, mWorkingRoot;
  protected boolean                  mVerbose;
  protected List<String>             mUsages;
  protected Config                   mConfig = ConfigFactory.load();
  //**************************************************************************
  // PRIVATE FIELDS
  //***-*********-*********-*********-*********-*********-*********-*********-
  protected String mAppName = "unTitled", mVersion, mBuild;

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
      if (args[i].length() == 0 || args[i].charAt(0) != '-') continue;

      // Is the argument a command-line option (starts with a '-') ?
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
      else if (args[i].equals("-w") || args[i].equals("-working"))
      {
        mWorkingRoot = args[i+1];
      }
      else if (args[i].equals("-o") || args[i].equals("-outfile"))
      {
        mOutfileRoot = args[i+1];
      }
      else if (args[i].equals("-cfg") || args[i].equals("-config"))
      {
        readConfig(args[i+1]);
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
    System.out.println("USAGE FOR " + mAppName + ": " + mBuild + " [options]");
    System.out.println("\n    COMMAND-LINE OPTIONS:");

    if (Tools.isSet(mUsages))
      for (String u : mUsages)
        System.out.println("    " + u);
  }
  abstract protected boolean doRun() throws Exception;
  abstract protected void    addProperty(String... vals);
  abstract public    String  getOutFile();

  public String getWorkingRoot() { return Strs.isSet(mWorkingRoot) ? mWorkingRoot+"/":""; }
  public String getLogFile()     { return getOutFile()+"_"+System.nanoTime()+".log"; }

  protected void close() throws IOException
  {
  }

  public String getAppName()         { return mAppName; }
  public void   setAppName(String s) { mAppName = s; }
  protected App addUsage(String tags, String usage)
  {
    if (mUsages==null) mUsages = new ArrayList<>();
    mUsages.add(tags + (Strs.isSet(usage) ? "\t:\t" + usage : ""));

    return this;
  }
  protected String option(String var, String usage, String[] args, int i, String... tags)
  {
    if (args[i].charAt(0)=='-' && Tools.isSet(args) && args.length>i+1 && Tools.isA(args[i], tags))
    {
      var=args[i+1];
      if (mUsages==null) mUsages = new ArrayList<>();
      mUsages.add(Strs.toString(tags, "\t") + (Strs.isSet(usage) ? "\t:\t" + usage : ""));
    }
    return var;
  }
  protected void addParamKey(String... names)
  {
    if (names!=null && names.length>1) sParamKeys.put(names[0], names[1]);
  }
  protected App readConfig(String cfgname)
  {
    try
    {
      System.out.println("Reading the configuration from " + cfgname);
      BufferedReader cfg = null;
      try
      {
        cfg = new BufferedReader(new InputStreamReader(new FileInputStream(cfgname)));
        while (cfg.ready())
        {
          String line = cfg.readLine().trim();
          // ignore the comments
          if (line.indexOf("//")==0 || line.indexOf("#")==0) continue;

          if (line.indexOf(":")>0) addParamKey(Strs.split(line, ':', true));
          else                     addProperty(Strs.split(line, '=', true));
        }
      }
      finally {
        if (cfg!=null) cfg.close();
      }
    }
    catch (IOException ie)
    {
      ie.printStackTrace();
    }
    return this;
  }
  protected void logn()         { log("\n"); }
  protected void logn(String s) { log(s+"\n"); }
  protected void log( String s)
  {
    try
    {
      FileWriter w = null;
      try
      {
        w = new FileWriter(getLogFile(), true);
        w.write(s);
      }
      finally { if (w!=null) w.close(); }
    }
    catch (Exception e) { System.out.print(s); }
  }

//  protected Config option(Config config, String usage, String[] args, int i, String... tags)
//  {
//    if (args[i].charAt(0)=='-' && Tools.isSet(args) && args.length>i+1 && Tools.isA(args[i], tags))
//    {
//      var=args[i+1];
//      if (mUsages==null) mUsages = new ArrayList<>();
//      mUsages.add(Strs.toString(tags, "\t") + (Strs.isSet(usage) ? "\t:\t" + usage : ""));
//    }
//    return var;
//  }
}