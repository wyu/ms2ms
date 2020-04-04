package org.ms2ms.algo.DIA;

import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.data.ms.SRMGroup;
import org.ms2ms.io.mzMLReader;

import java.io.FileWriter;
import java.io.IOException;

public class DIA_utils
{
  public static MultiTreeTable<Float, Float, SRMGroup> runDIA(
      String root, String tr, String tag, Tolerance tol, Float rt_tol, String run) throws IOException
  {
    MultiTreeTable<Float, Float, SRMGroup> groups =  SRMGroup.readTransitions(root+tr);

    groups = mzMLReader.extractTransitionXICs(root, run+".mzML", tol, rt_tol, groups);

    FileWriter xic = new FileWriter(root+run+tag+".xic"), ftr = new FileWriter(root+run+tag+".feature");

    SRMGroup.headerXIC(xic); SRMGroup.headerFeatures(ftr);
    for (SRMGroup grp : groups.values())
    {
      grp.composite().centroid().scoreSimillarity();
      grp.printXIC(xic).printFeatures(ftr);
    }
    xic.close(); ftr.close();

    return groups;
  }
}
