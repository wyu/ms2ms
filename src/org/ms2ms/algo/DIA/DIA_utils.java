package org.ms2ms.algo.DIA;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.data.ms.SRMGroup;
import org.ms2ms.io.mzMLReader;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

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
      grp.composite().centroid(5f, 0.5f).scoreSimillarity().calcFeatureExclusivity(0.25f, 3);
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
      grp.impute(lc_width*0.1f).composite().centroid(5f, lc_width).scoreSimillarity().calcFeatureExclusivity(lc_width, 3);
      grp.printXIC(xic).printFeatures(ftr);
    }

    // gather the protei/peptides
    Multimap<String, String> proteins = HashMultimap.create();
    if (Tools.isSet(peptides))
      for (String s : peptides) proteins.put("unTitled", s);

    for (SRMGroup grp : groups.values())
      if (grp.getProteinId()!=null) proteins.put(grp.getProteinId(), grp.getSequence());

    if (Tools.isSet(proteins))
      for (String pid : proteins.keySet())
      {
        SRMGroup protein = SRMGroup.buildProteinProfile(groups.values(), pid, true, Strs.toStringArray(proteins.get(pid)));

        protein.centroid(5f, lc_width*60f).calcFeatureExclusivity(lc_width*60f, 3); // in secs
        protein.printXIC(xic).printFeatures(ftr);
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
