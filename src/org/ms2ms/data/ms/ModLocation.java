package org.ms2ms.data.ms;

import org.ms2ms.utils.Tools;

/**
 * Created by yuw on 8/15/16.
 */
public class ModLocation implements Cloneable
{
  public int    locations;
  public double mods=0, occupancy=0;
  public NLoss nloss=null;

  public ModLocation(int loc, double mm)
  {
    locations=loc; mods=mm;
  }

  public ModLocation setOccupancy(double s) { occupancy=s; return this; }
  public ModLocation setLocation(    int s) { locations=s; return this; }
  public ModLocation setNLoss(     NLoss s) { nloss    =s; return this; }

  @Override
  public ModLocation clone()
  {
    ModLocation cloned = new ModLocation(locations, mods);
    if (nloss!=null) cloned.setNLoss(nloss.clone());
    return cloned;
  }
  @Override
  public String toString()
  {
    return locations + "@" + Tools.d2s(mods, 3);
  }
}
