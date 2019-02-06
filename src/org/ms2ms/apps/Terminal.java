package org.ms2ms.apps;

import java.io.FileWriter;
import java.io.PrintStream;

public class Terminal
{
  protected PrintStream mConsole = new PrintStream(System.out, true);
  protected String mLogFile = null, mJobID="";
  protected boolean mIsTestOnly=false;

  public Terminal setLogfile(String s) { mLogFile=s; return this; }

  public void logn()           { log("\n", false); }
  public void printn()         { log("\n", true); }
  public void log(   String s) { log(s, false); }
  public void logn(  String s) { log(s+"\n", false); }
  public void logn(  StringBuffer s) { if (s!=null) log(s.toString()+"\n", false); }
  public void printn(StringBuffer s) { if (s!=null) log(s.toString()+"\n", true); }
  public void printn(String s) { log(s+"\n", true); }
  public void print(String s) { log(s, true); }

  public void log( String s, boolean verbose)
  {
    if ("NIL".equals(mLogFile))
    {
      if (verbose) mConsole.print(s);
      return;
    }

    try
    {
      FileWriter w = null;
      try
      {
        w = new FileWriter(mLogFile, true);
        w.write(s);
        if (verbose) mConsole.print(s);
      }
      finally { if (w!=null) w.close(); }
    }
    catch (Exception e) { mConsole.print(s); }
  }

  public boolean isTestOnly()          { return mIsTestOnly; }
  public void    isTestOnly(boolean s) { mIsTestOnly=s; }
}
