package org.ms2ms.algo.DIA;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.data.ms.SRMGroup;
import org.ms2ms.io.mzMLReader;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

public class DIA_utils
{
  public static MultiTreeTable<Float, Float, SRMGroup> runDIA(
      MultiTreeTable<Float, Float, SRMGroup> groups,
      String root, String tag, Tolerance tol, float rt_tol, String run) throws IOException
  {
//    MultiTreeTable<Float, Float, SRMGroup> groups =  SRMGroup.readTransitions(root+tr);

    groups = mzMLReader.extractTransitionXICs(root, run+".mzML", tol, rt_tol, groups);

    FileWriter xic = new FileWriter(root+run+tag+".xic"), ftr = new FileWriter(root+run+tag+".feature");

    SRMGroup.headerXIC(xic); SRMGroup.headerFeatures(ftr);
    for (SRMGroup grp : groups.values())
    {
      grp.composite().centroid(5f, 0.5f).scoreSimillarity().calcFeatureExclusivity(0.25f);
      grp.printXIC(xic).printFeatures(ftr);
    }
    xic.close(); ftr.close();

    return groups;
  }
  public static ListMultimap<Integer, Float> buildFragmentBank(Collection<SRMGroup> groups, float ratio)
  {
    ListMultimap<Integer, Float> bank = ArrayListMultimap.create();
    for (SRMGroup grp : groups)
      for (Float frag : grp.getSRMs().keySet()) bank.put((int )Math.round(frag*0.01), frag);

    return bank;
  }
}
