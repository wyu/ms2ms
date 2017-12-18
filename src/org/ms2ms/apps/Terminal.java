package org.ms2ms.apps;

import java.io.FileWriter;

public class Terminal
{
  protected String mLogFile = null, mJobID="";

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
      if (verbose) System.out.print(s);
      return;
    }

    try
    {
      FileWriter w = null;
      try
      {
        w = new FileWriter(mLogFile, true);
        w.write(s);
        if (verbose) System.out.print(s);
      }
      finally { if (w!=null) w.close(); }
    }
    catch (Exception e) { System.out.print(s); }
  }
}
