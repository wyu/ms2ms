package org.ms2ms.test.io;

import com.google.common.collect.Multimap;
import de.mpc.pia.intermediate.compiler.PIACompiler;
import org.apache.commons.io.filefilter.WildcardFileFilter;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import org.junit.Test;
import org.ms2ms.algo.PSMs;
import org.ms2ms.data.ms.Engine;
import org.ms2ms.io.MsfReader;
import org.ms2ms.io.PsmReaders;
import org.ms2ms.io.PsmWriters;
import org.ms2ms.test.TestAbstract;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;

import java.util.List;

/**
 * ** Copyright 2014-2015 ms2ms.org
 * <p/>
 * Description:
 * <p/>
 * Author: wyu
 * Date:   2/10/15
 */
public class IOsTest extends TestAbstract
{
//  @Test
//  public void msf() throws Exception
//  {
//    PIACompiler pia = new PIACompiler();
//    String msfile = "/Users/yuw/Documents/Data/PIIQ/151124_DBH_Mosser_TMT1_MS2_All_V1.msf";
//    MsfReader.getDataFromThermoMSFFile("test", "", pia);
//  }
  @Test
  public void mzID2Novor() throws Exception
  {
    String root = "/Users/yuw/Documents/Apps/Joslin/PSMs/";
    Multimap<SpectrumIdentifier, PeptideMatch> matches = PsmReaders.readAmanda(root+"default/Amanda_mgf_output_2.csv", 1);
    PsmWriters.writeNovor(root + "F01.novor.csv", matches, Engine.AMANDA.getCanonicalScore(), "Mouse_Plasma_LIRKO2_01_27Aug12_Lynx_12-06-05");
  }
  @Test
  public void StrsSplit() throws Exception
  {
    String p1 = "IETLMRNLM[15.9949]PWRK",
           p2 = "+229.163HMK+229.163K+229.163HAK+229.163K+229.163MK+229.163K+229.163QMK+229.163K+229.163";
    String r1 = "(^.*?\\[|\\]\\s*$)", r2 = "\\]\\s*,\\s*\\[";

    PeptideMatch m1 = PSMs.fromNumModSequence(p1),
                 m2 = PSMs.fromNumModSequence(p2);
    List<String> item01 = Strs.splits(p1, "[+-.\\d]+"),
    item02 = Strs.splits(p1, r1),
    item03 = Strs.splits(p1, r2);

    System.out.println(item01.size());
  }
  @Test
  public void recursiveListing()
  {
    IOs.listFiles("/tmp", new WildcardFileFilter("*.tmp"));
  }
}
