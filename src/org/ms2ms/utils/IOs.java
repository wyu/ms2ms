package org.ms2ms.utils;

import com.google.common.collect.Table;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 7/13/14
 * Time: 3:50 PM
 * To change this template use File | Settings | File Templates.
 */
public class IOs
{
  public static boolean exists(String s)
  {
    if (!Tools.isSet(s)) return false;
    return new File(s).exists();
  }
  public static void write(String f, String data)
  {
    FileWriter w = null;
    try
    {
      try
      {
        w = new FileWriter(f);
        w.write(data);
      }
      catch (IOException ie)
      {
        System.out.println("Failed to write the output to " + f);
      }
      finally
      {
        if (w!=null) w.close();
      }
    }
    catch (IOException ie) {}
  }
}
