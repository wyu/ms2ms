package org.ms2ms.data.ms;

public class NLoss implements Cloneable
{
    private String mName;
    private double nloss=0, ratio=0;

    public NLoss()         { mName="unTitled"; }
    public NLoss(String s) { mName=s; }
    public NLoss(String s, double nl, double r) { nloss=nl; ratio=r; mName=s; }

    public double getNLoss() { return nloss; }
    public double getRatio() { return ratio; }
    public String getName()  { return mName; }

    @Override
    public NLoss clone()
    {
      return new NLoss(mName, nloss, ratio);
    }
}
