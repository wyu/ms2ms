package org.ms2ms.test.io;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.apache.commons.io.filefilter.WildcardFileFilter;
import org.expasy.mzjava.proteomics.io.mol.FastaProteinReader;
import org.expasy.mzjava.proteomics.mol.Protein;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import org.junit.Test;
import org.ms2ms.algo.PSMs;
import org.ms2ms.data.ms.Engine;
import org.ms2ms.io.CustomFastaProteinReader;
import org.ms2ms.io.MsIO;
import org.ms2ms.io.PsmReaders;
import org.ms2ms.io.PsmWriters;
import org.ms2ms.test.TestAbstract;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.TabFile;
import org.ms2ms.utils.Tools;

import java.io.*;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

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
  public void lineByline() throws Exception
  {
    String root = "/Users/kfvf960/Apps/pipeline/Chorus/sequences/Uniprot/";
    // go thro the file line by line to correct syntactical miskate
    BufferedReader R = IOs.newBufferedReader(root+"uniprot_trembl.xml.gz", "UTF-8");

    OutputStream fileStream = new GZIPOutputStream(new FileOutputStream(root+"uniprot_trembl_fixed.xml.gz"));
    BufferedWriter W = new BufferedWriter(new OutputStreamWriter(fileStream, "UTF-8"));

    System.out.println("Fixing the lines...");
    long rows=0;
    while (R.ready())
    {
      String line = R.readLine();
      //</entry><entry dataset="Swiss-Prot" created="2009-06-16" modified="2017-09-27" version="18" xmlns="http://uniprot.org/uniprot">
      if (line.indexOf("</entry><entry")==0)
        line = line.replaceAll("</entry><entry","</entry>\n<entry");
      else if (line.indexOf("</entry></uniprot>")==0)
        line = line.replaceAll("</entry></uniprot>","</entry>\n</uniprot>");

      W.write(line+"\n");

      if (++rows % 100000==0) System.out.print(".");
      if (rows % 10000000==0) System.out.println(rows);
    }
    R.close();
    W.close();
    fileStream.close();
  }
  @Test
  public void idmapping() throws Exception
  {
    // http://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/
    Multimap<String, String> gene_cho = HashMultimap.create();

    String root = "/Users/kfvf960/Apps/pipeline/KB/idmapping/";
    TabFile cho2mouse = new TabFile(root+"CHO_Mouse_Rat_Ortho_20180207.txt", TabFile.tabb);
    while (cho2mouse.hasNext())
    {
      // Gene stable ID	Transcript stable ID	Gene name	Golden Hamster gene stable ID	Golden Hamster gene name	Mouse gene stable ID	Mouse gene name	Mouse protein or transcript stable ID	Golden Hamster protein or transcript stable ID	Rat gene name	Rat gene stable ID	Rat protein or transcript stable ID
      // ENSCGRG00001000006	ENSCGRT00001000006	NADH1	ENSMAUG00000000006	ND1	ENSMUSG00000064341	mt-Nd1	ENSMUSP00000080991	ENSMAUP00000000001	Mt-nd1	ENSRNOG00000030644	ENSRNOP00000049172

      if (Strs.isSet(cho2mouse.get("Gene name")) && Strs.isSet(cho2mouse.get("Golden Hamster gene name")))
      {
        gene_cho.put(cho2mouse.get("Gene name"), cho2mouse.get("Golden Hamster gene name"));
      }
    }
    cho2mouse.close();
    FileWriter w = new FileWriter(root+"CHO_Mouse_Rat_20180207.txt");
    w.write("CHOK1\tHuman\n");
    for (String gene : gene_cho.keySet())
      for (String cho : gene_cho.get(gene))
        w.write(cho+"\t"+gene+"\n");
    w.close();
  }
  @Test
  public void asUniprotFASTA() throws Exception
  {
    String                db = "/Users/kfvf960/Apps/pipeline/Chorus/sequences/", tag="rnaseq", species="Homo sapiens";
    FastaProteinReader fasta = new CustomFastaProteinReader(new File(db+"merged.pep.fasta"));
    FileWriter             w = new FileWriter(db+"HeLa_pep_UP.fasta");
    while (fasta.hasNext())
    {
//      >MT.A2M.ENST00000318602.missense.749A/S
//      LVHVEEPHTETVRKYFPETWIWDLVVVNSSGVAEVGVTVPDTITEWKAGAFCLSEDAGL
      Protein p = fasta.next();

//      >sp|A0A075B6I1|LV460_HUMAN Immunoglobulin lambda variable 4-60 OS=Homo sapiens GN=IGLV4-60 PE=3 SV=1
//      MAWTPLLLLFPLLLHCTGSLSQPVLTQSSSASASLGSSVKLTCTLSSGHSSYIIAWHQQQ
//      PGKAPRYLMKLEGSGSYNKGSGVPDRFSGSSSGADRYLTISNLQFEDEADYYCETWDSNT
//      >tr|A0A075B6I3|A0A075B6I3_HUMAN Immunoglobulin lambda variable 11-55 (non-functional) OS=Homo sapiens GN=IGLV11-55 PE=4 SV=2
//      MALTPLLLLLLSHCTGSLSRPVLTQPPSLSASPGATARLPCTLSSDLSVGGKNMFWYQQK
//      LGSSPRLFLYHYSDSDKQLGPGVPSRVSGSKETSSNTAFLLISGLQPEDEADYYCQVYES
//      SANHSETDEEVGQKPRF

      String[] tokens = Strs.split(p.getAccessionId(), '.', true);
      String ac = p.getAccessionId(), name = p.getAccessionId(), misc =" OS="+species+" GN="+tokens[1];
      // - Warning - accession [MT.ABCC1.ENST00000346370.FS.1422-1426GGACTGCACCGT/G] is longer than 50 characters

      if (ac.length()>49)
      {
        ac = "";
        for (int k=tokens.length-1; k>=0; k--)
        {
          if (ac.length()+tokens[k].length()>48) break;
          ac = Strs.extend(tokens[k],ac, ".");
        }
      }

      w.write(">"+tag+"|"+ac+" "+name+misc+"\n");
      for (int i=0; i<p.size(); i+=60)
      {
        w.write(p.toSymbolString().substring(i, Math.min(i+60, p.size()))+"\n");
        if (i+60>p.size()) break;
      }
    }
    fasta.close(); w.close();
  }
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
    IOs.listFiles("/tmp", new WildcardFileFilter("*.tmp"), 1);

    List<String> files = MsIO.listFiles("/Users/kfvf960/Apps/pipeline/Chorus/data/*Scan*.txt");
    System.out.println(files.size());
  }
}
