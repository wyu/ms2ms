package org.ms2ms.algo.DIA;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.data.ms.SRM;
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

    groups = mzMLReader.extractTransitionXICs(root, run+".mzML", tol, rt_tol, groups, 0.5f);

    FileWriter xic = new FileWriter(root+run+tag+".xic"), ftr = new FileWriter(root+run+tag+".feature");

    SRMGroup.headerXIC(xic); SRMGroup.headerFeatures(ftr);
    for (SRMGroup grp : groups.values())
    {
      grp.composite().centroid(5f, 0.5f, false).scoreSimillarity().calcFeatureExclusivity(0.25f, 3);
      grp.printXIC(xic).printFeatures(ftr);
    }
    xic.close(); ftr.close();

    return groups;
  }
  public static MultiTreeTable<Float, Float, SRMGroup> runPRM(
      MultiTreeTable<Float, Float, SRMGroup> groups,
      String run, String root, String tag, Tolerance tol, float rt_tol, float lc_width, String... peptides) throws IOException
  {
    System.out.println(run);
//    MultiTreeTable<Float, Float, SRMGroup> groups =  SRMGroup.readTransitions(root+tr);

    groups = mzMLReader.extractTransitionXICs(root, run+".mzML", tol, rt_tol, groups, 0f);

    FileWriter xic = new FileWriter(root+run+tag+".xic"), ftr = new FileWriter(root+run+tag+".feature");

    SRMGroup.headerXIC(xic); SRMGroup.headerFeatures(ftr);
    for (SRMGroup grp : groups.values())
    {
      grp.impute(lc_width*0.1f).composite().centroid(5f, lc_width, false).scoreSimillarity().calcFeatureExclusivity(0.25f, 3);
      grp.printXIC(xic).printFeatures(ftr);
    }

    SRM profile = SRMGroup.composite(groups.values(), peptides);
    SRMGroup protein = new SRMGroup("Protein Summary");

    protein.getSRMs().put(0f, profile);

    protein.centroid(5f, lc_width, true);
    protein.printXIC(xic).printFeatures(ftr);

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
