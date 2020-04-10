package org.ms2ms.test.ms;

import com.google.common.collect.ListMultimap;
import org.expasy.mzjava.core.ms.Tolerance;
import org.junit.Test;
import org.ms2ms.algo.DIA.DIA_utils;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.data.ms.OffsetPpmTolerance;
import org.ms2ms.data.ms.SRMGroup;
import org.ms2ms.io.mzMLReader;
import org.ms2ms.r.Dataframe;
import org.ms2ms.test.TestAbstract;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

public class DIATest extends TestAbstract
{
  //  String root = "/Users/hliu/Desktop/App/2014/data/mzXML-centroid/";
  String root = "/Users/kfvf960/OneDrive - AZCollaboration/TechDev/DIA/PoC_Alignment/XIC/";

  // TODO: Need to modify or extend MzxmlReader to readSpectrumIdentifier only selected msLevel or RT range, etc
  // peak processing takes lots of time!

  @Test
  public void runLandmarks() throws Exception
  {
    MultiTreeTable<Float, Float, SRMGroup> groups =  SRMGroup.readTransitions(root+"L039_SERD_Project_4Frag_007_tr.tsv");

    DIA_utils.runDIA(groups, root, "_landmarks_30sec", new OffsetPpmTolerance(15d), 5f, "M1A1732_20191211_SERD_AA178");
  }
  @Test
  public void runLibrary() throws Exception
  {
    MultiTreeTable<Float, Float, SRMGroup> groups =  SRMGroup.readTransitions(root+"L039_SERD_Project_4Frag_007_lib_tr.tsv");

//    DIA_utils.runDIA(groups, root, "_lib", new OffsetPpmTolerance(15d), 5f, "M1A1732_20191211_SERD_AA178");

    // turn the group into decoys
    ListMultimap<Integer, Float> bank = DIA_utils.buildFragmentBank(groups.values(), 0.01f);
    for (SRMGroup group : groups.values()) group.mutate(bank, new Random(System.nanoTime()));

    DIA_utils.runDIA(groups, root, "_lib_mutated", new OffsetPpmTolerance(15d), 5f, "M1A1732_20191211_SERD_AA178");
  }

}
