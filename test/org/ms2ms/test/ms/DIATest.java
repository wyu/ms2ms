package org.ms2ms.test.ms;

import org.junit.Test;
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
  public void readTransitionList() throws Exception
  {
    MultiTreeTable<Float, Float, SRMGroup> groups =  SRMGroup.readTransitions(root+"L039_SERD_Project_4Frag_007_tr.tsv");

    OffsetPpmTolerance tol = new OffsetPpmTolerance(15d);
    groups = mzMLReader.extractTransitionXICs(root, "M1A1732_20191211_SERD_AA178.mzML", tol, 5f, groups);

    FileWriter xic = new FileWriter(root+"M1A1732_20191211_SERD_AA178.xic"),
               ftr = new FileWriter(root+"M1A1732_20191211_SERD_AA178.feature");

    SRMGroup.headerXIC(xic); SRMGroup.headerFeatures(ftr);
    for (SRMGroup grp : groups.values())
    {
      grp.composite().printXIC(xic);
      grp.centroid().printFeatures(ftr);
    }
    xic.close(); ftr.close();
  }

}
