package org.ms2ms.test.ms;

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

public class DIATest extends TestAbstract
{
  //  String root = "/Users/hliu/Desktop/App/2014/data/mzXML-centroid/";
  String root = "/Users/kfvf960/OneDrive - AZCollaboration/TechDev/DIA/PoC_Alignment/XIC/";

  // TODO: Need to modify or extend MzxmlReader to readSpectrumIdentifier only selected msLevel or RT range, etc
  // peak processing takes lots of time!

  @Test
  public void runLandmarks() throws Exception
  {
    DIA_utils.runDIA(root, "L039_SERD_Project_4Frag_007_tr.tsv", "_landmarks_30sec", new OffsetPpmTolerance(15d), 5f, "M1A1732_20191211_SERD_AA178");
  }
  @Test
  public void runLibrary() throws Exception
  {
    DIA_utils.runDIA(root, "L039_SERD_Project_4Frag_007_lib_tr.tsv", "_lib", new OffsetPpmTolerance(15d), 5f, "M1A1732_20191211_SERD_AA178");
  }

}
