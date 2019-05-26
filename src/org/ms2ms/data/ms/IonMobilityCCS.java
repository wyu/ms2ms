package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.spectrum.RetentionTime;

public class IonMobilityCCS implements RetentionTime
{
  private final double ccs;

  public IonMobilityCCS(double s) { this.ccs = s; }
  public IonMobilityCCS(IonMobilityCCS src) {
    this.ccs = src.ccs;
  }

  public double getCCS() { return ccs; }

  @Override
  public double getTime() { return ccs; }

  @Override
  public double getMinRetentionTime() { return ccs; }

  @Override
  public double getMaxRetentionTime() { return ccs; }

  @Override
  public RetentionTime copy() { return new IonMobilityCCS(this); }

  @Override
  public boolean contains(RetentionTime src) { return this.getTime() == src.getTime(); }
}
