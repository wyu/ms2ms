package org.ms2ms.test.ms;

import com.google.common.collect.Multimap;
import org.junit.Test;
import org.ms2ms.algo.Isobarics;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.data.ms.PeptideFeature;
import org.ms2ms.r.Dataframe;
import org.ms2ms.test.TestAbstract;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.TabFile;

import java.io.FileWriter;

/**
 * ** Copyright 2014-2015 ms2ms.org
 * <p/>
 * Description:
 * <p/>
 * Author: wyu
 * Date:   3/26/15
 */
public class IsobaricsTest extends TestAbstract
{

  @Test
  public void QQPlot() throws Exception
  {
    MultiTreeTable<String, String, PeptideFeature> peptides = newTMT();

    Multimap<String, String> grps = Strs.newMultimap('=',"Ctrl=S1;S2;S3;S4;S5;S6;F1;F2;F3;F4","COPD=CP1;CP2;CP3;CP4;CP5;CP6;CP7;CP8;CP9;CP10");
    Dataframe data = Isobarics.fromGeneExpt(peptides, grps);

    FileWriter file = new FileWriter("/home/wyu/Projects/reports/data.frm");
    file.write(data.display().toString());
    file.close();
}
  public static MultiTreeTable<String, String, PeptideFeature> newTMT() throws Exception
  {
    MultiTreeTable<String, String, PeptideFeature> gene_seq_peptide = Isobarics.Scan2GeneExpt("set1", null,
        "GeneName", "Peptide Sequence", new TabFile("/home/wyu/Projects/ms2ms/data/TMT_SET1.txt", TabFile.tabb),
        Strs.newMap('=', "c127 (CP1)=CP1", "c128 (CP2)=CP2", "c129 (CP6)=CP6", "c130 (S1)=S1", "c131 (S2)=S2"), "c126 (CC1)");

    System.out.println(gene_seq_peptide.getData().keySet().size() + "/" + gene_seq_peptide.size());
    gene_seq_peptide = Isobarics.Scan2GeneExpt("set2", gene_seq_peptide,
        "GeneName", "Peptide Sequence", new TabFile("/home/wyu/Projects/ms2ms/data/TMT_SET2.txt", TabFile.tabb),
        Strs.newMap('=', "c127 (CP8)=CP8","c128 (CP9)=CP9","c129 (CP10)=CP10","c130 (S3)=S3","c131 (S4)=S4"), "c126 (CC2)");

    System.out.println(gene_seq_peptide.getData().keySet().size() + "/" + gene_seq_peptide.size());
    gene_seq_peptide = Isobarics.Scan2GeneExpt("set3", gene_seq_peptide,
        "GeneName", "Peptide Sequence", new TabFile("/home/wyu/Projects/ms2ms/data/TMT_SET3.txt", TabFile.tabb),
        Strs.newMap('=', "c127 (CP5)=CP5","c128 (CP7)=CP7","c129 (F1)=F1","c130 (F4)=F4","c131 (S5)=S5"), "c126 (CC3)");

    System.out.println(gene_seq_peptide.getData().keySet().size() + "/" + gene_seq_peptide.size());
    gene_seq_peptide = Isobarics.Scan2GeneExpt("set4", gene_seq_peptide,
        "GeneName", "Peptide Sequence", new TabFile("/home/wyu/Projects/ms2ms/data/TMT_SET4.txt", TabFile.tabb),
        Strs.newMap('=', "c127 (S6)=S6","c128 (CP4)=CP4","c129 (F2)=F2","c130 (F3)=F3","c131 (CP3)=CP3"), "c126 (CC4)");

    System.out.println(gene_seq_peptide.getData().keySet().size() + "/" + gene_seq_peptide.size());

    return gene_seq_peptide;
  }
  public static MultiTreeTable<String, String, PeptideFeature> newTMT0() throws Exception
  {
    MultiTreeTable<String, String, PeptideFeature> gene_seq_peptide = Isobarics.Scan2GeneExpt("unTitled", null,
        "GeneName", "Peptide Sequence", new TabFile("/home/wyu/Projects/ms2ms/data/TMT_SET1_tops.csv", TabFile.tabb),
        Strs.newMap('=', "c127 (CP1)=CP1", "c128 (CP2)=CP2", "c129 (CP6)=CP6", "c130 (S1)=S1", "c131 (S2)=S2"), "c126 (CC1)");

    return gene_seq_peptide;
  }
}
