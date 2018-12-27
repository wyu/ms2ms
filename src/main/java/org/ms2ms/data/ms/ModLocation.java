package org.ms2ms.data.ms;

import org.ms2ms.utils.Tools;

/**
 * Created by yuw on 8/15/16.
 */
public class ModLocation
{
  public int    locations;
  public double mods;

  public ModLocation(int loc, double mm)
  {
    locations=loc; mods=mm;
  }

  @Override
  public String toString()
  {
    return locations + "@" + Tools.d2s(mods, 3);
  }
}
