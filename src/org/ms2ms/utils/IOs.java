package org.ms2ms.utils;

import com.google.common.collect.Table;

import java.io.File;

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
}