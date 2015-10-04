package org.ms2ms.algo;

import com.google.common.collect.Multimap;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.data.ms.PeptideFeature;
import org.ms2ms.r.Dataframe;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.TabFile;
import org.ms2ms.utils.Tools;

import java.io.IOException;
import java.util.Collection;
import java.util.Map;

/** For iTraq or TMT labeling experiment
 * ** Copyright 2014-2015 ms2ms.org
 * <p/>
 * Description:
 * <p/>
 * Author: wyu
 * Date:   3/26/15
 */
public class Isobarics
{
  public static final String TITLE    = "title";
  public static final String GENE     = "gene";
  public static final String PEPTIDE  = "peptide";
  public static final String CHARGE   = "z";
  public static final String GROUP    = "group";
  public static final String SAMPLE   = "sample";
  public static final String ABUND    = "abund";
  public static final String MEAN     = "avg";
  public static final String RAW_MEAN = "raw.avg";
  public static final String BASIS    = "basis";

  /** Transform the peptide table into a gene-by-expt table. The row is stored as
   *  a col-map in the cell
   *
   *
   *
   * @param title
   * @param set is the tabular file with the data
   * @param col_expts is an array of the 'channels' corresponding to the isobaric labels
   * @return a table with the transformed data
   *
   * @throws IOException
   */
  public static MultiTreeTable<String, String, PeptideFeature>
      Scan2GeneExpt(String title, MultiTreeTable<String, String, PeptideFeature> gene_seq_val,
                    String col_gene, String col_seq,
                    TabFile set, Map<String, String> col_expts, String basis) throws IOException
  {
//    Accession Number	Protein Redundancy	Protein Description	Protein Masses	Peptide Sequence	Variable Modifications	Protein-Relative Modifications	Experimental mz	Charge	Predicted mr	Delta	Peptide Score	Start Position	End Position	Preceding Residue	Following Residue	Missed Cleavages	GeneName	c126 (CC1)	c127 (CP1)	c128 (CP2)	c129 (CP6)	c130 (S1)	c131 (S2)
//    gi|254910983|ref|NP_775782.2|	1	"amyotrophic lateral sclerosis 2 (juvenile) chromosome region, candidate 13 [Homo sapiens]"	69146	RNGSPTPAGSLGGGAVATAGGPGSR	S4: Phospho	S10: Phospho	821.07914	3	2460.207367	0.008225	126.29	7	31	R	L	1	['FAM117B']	1.62E+03	2.91E+03	1.66E+03	0.00E+00	2.28E+03	1.86E+03
//    gi|23821044|ref|NP_667339.1|	1	zinc finger/RING finger 2 [Homo sapiens]	26217	AYSGSDLPSSSSGGANGTAGGGGGAR	S3: Phospho	S19: Phospho	836.37224	3	2506.092453	0.002439	119.93	17	42	R	A	0	['ZNRF2']	1.62E+03	1.22E+03	2.11E+03	0.00E+00	1.68E+03	2.29E+03
    if (gene_seq_val==null) gene_seq_val = MultiTreeTable.create();
    while (set.hasNext())
    {
      String gene = Strs.isSet(set.get(col_gene))?Strs.split(set.get(col_gene),'\'')[1]:"DUMMY";
      Double base = set.getDouble(basis);
      if (base!=null && base!=0)
        // create a new peptide feature for the row
        gene_seq_val.put(gene, set.get(col_seq),
            PeptideFeature.fromMultiplierZ(title, set.getMappedRow(), col_expts, base));
    }
    return gene_seq_val;
  }
  public static Dataframe fromGeneExpt(MultiTreeTable<String, String, PeptideFeature> set, Multimap<String, String> grps)
  {
    if (!Tools.isSet(set)) return null;

    Dataframe data = new Dataframe(); int rows=0;
    // for each gene, peptide combination
    for (String gene : set.getData( ).keySet())
    {
      if (gene.indexOf("unresolved")==0) continue;
      for (String seq : set.row(gene).keySet())
      {
        Collection<PeptideFeature> peptides = set.get(gene, seq);
        // get the group average
        for (String grp : grps.keySet())
        {
          double avg = 0d, raw=0, n=0;
          for (PeptideFeature peptide : peptides)
            // going thro each experimental group
            for (String expt : grps.get(grp))
              if (peptide.getAbundance(expt)!=null)
              {
                avg+=peptide.getAbundance(expt);
                raw+=peptide.getAbundance(expt)*peptide.getAbundance(BASIS);
                n++;
              }
          avg/=n; raw/=n;

          for (String expt : grps.get(grp))
            for (PeptideFeature peptide : peptides)
              if (peptide.getAbundance(expt)!=null)
              {
                data.put(++rows,GENE,gene).put(PEPTIDE,seq).put(GROUP,grp).put(SAMPLE,expt).put(MEAN,avg).put(TITLE,peptide.getTitle());
                data.put(ABUND, peptide.getAbundance(expt)).put(BASIS, peptide.getAbundance(BASIS)).put(RAW_MEAN,raw);
              }
        }
      }
    }

    return data.init();
  }
}
