package org.ms2ms.data.ms;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by yuw on 8/20/16.
 */
public class ModSetting
{
  public double[]          masses;
  public List<ModLocation> mods;

  public ModSetting(int size, int loc, double increment)
  {
    masses = new double[size];
    mods = new ArrayList<>(); mods.add(new ModLocation(loc, increment));
  }
  public ModSetting(double[] m, List<ModLocation> mod)
  {
    masses=m; mods=mod;
  }
}
