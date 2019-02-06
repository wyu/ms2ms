package org.ms2ms.test.io;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.expasy.mzjava.proteomics.io.mol.FastaProteinReader;
import org.expasy.mzjava.proteomics.mol.Protein;
import org.junit.Test;
import org.ms2ms.io.CustomFastaProteinReader;
import org.ms2ms.test.TestAbstract;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

public class ProteinSequenceTest extends TestAbstract
{
  public final String NOTMAPPED = "NotMapped";

  @Test
  public void FASTA2Tab() throws Exception
  {
    String root = "/Users/kfvf960/Apps/pipeline/Chorus/sequences/";
    Multimap<String, Protein> seqs = readProteinSequences(
        root+"48ProteinUPS1.fasta", root+"48ProteinUPS1_Bonus.fasta", root+"contaminants.fasta",
        root+"craps.fasta", root+"vectors.fasta", root+"proteases.fasta");

    writeProteinsAsTab(root+"std48_expected.tab", seqs.values());
  }
  @Test
  public void consolidate() throws Exception
  {
    String root = "/Users/kfvf960/Apps/pipeline/Chorus/sequences/";
    Multimap<String, Protein> originals = readProteinSequences(root+"contaminants.fasta", root+"craps.fasta"),
                              uniprots  = readProteinSequences(root+"contaminants_uniprot_isoforms.fasta", root+"craps_uniprot.fasta");

    // figure out the mapping
    Multimap<String, Protein> mapped = map2uniprot(originals, uniprots);
    if (mapped.containsKey(NOTMAPPED))
    {
      writeProteinsAsFASTA(root+"notmapped.fasta", mapped.get(NOTMAPPED));
      mapped.removeAll(NOTMAPPED);
    }
    if (Tools.isSet(mapped))
      writeProteinsAsFASTA(root+"mapped.fasta", mapped.values());
  }
  @Test
  public void FASTA2SeqLines() throws Exception
  {
//    StringBuffer seqlines = FASTA2Lines("/Users/kfvf960/Apps/pipeline/Chorus/sequences/contaminants.fasta", 100);
    StringBuffer seqlines = FASTA2Lines("/Users/kfvf960/Apps/pipeline/Chorus/sequences/craps.fasta", 100);
    IOs.save("/tmp/seqs.txt", seqlines);
  }
  private void writeProteinsAsFASTA(String outfile, Collection<Protein> proteins) throws IOException
  {
    FileWriter w = new FileWriter(outfile);
    for (Protein p : proteins)
    {
      w.write(">"+p.getAccessionId()+"\n");
      w.write(p.toSymbolString()+"\n");
    }
    w.close();
  }
  private void writeProteinsAsTab(String outfile, Collection<Protein> proteins) throws IOException
  {
    FileWriter w = new FileWriter(outfile);

    w.write("ID\tSequence\n");
    for (Protein p : proteins)
    {
      w.write(p.getAccessionId()+"\t");
      w.write(p.toSymbolString()+"\n");
    }
    w.close();
  }
  private Multimap<String, Protein> map2uniprot(Multimap<String, Protein> originals, Multimap<String, Protein> uniprots)
  {
    Multimap<String, Protein> mapped = HashMultimap.create();
    // figure out the mapping
    for (String original : originals.keySet())
    {
      String found = Strs.getEqualOrEmbedded(original, uniprots.keySet());
      Protein sw=null, tr=null;
      if (found!=null)
      {
        System.out.println("Looking for the best mapping...");
        for (Protein p1 : uniprots.get(original))
        {
          if (sw==null && p1.getAccessionId().indexOf("sp")==0) sw=p1;
          if (tr==null && p1.getAccessionId().indexOf("tr")==0) tr=p1;
        }
      }
      if      (sw!=null) mapped.put(sw.toSymbolString(), sw);
      else if (tr!=null) mapped.put(tr.toSymbolString(), tr);
      else
      {
        System.out.print("Not mapped! Sent to manual mapping");
        for (Protein p0 : originals.get(original)) System.out.print(p0.getAccessionId()+"; ");
        System.out.println();
        mapped.put(NOTMAPPED, Tools.front(originals.get(original)));
      }
    }

    return mapped;
  }
  private StringBuffer FASTA2Lines(String db, int bundle) throws IOException
  {
    FastaProteinReader fasta = new CustomFastaProteinReader(new File(db));
    StringBuffer         buf = new StringBuffer();

    int rows=0;
    while (fasta.hasNext())
    {
      Protein p = fasta.next();
      buf.append(p.toSymbolString()+"\n");
      if (++rows % bundle==0) buf.append("\n");
    }
    fasta.close();

    return buf;
  }
  private Multimap<String, Protein> readProteinSequences(String... files) throws IOException
  {
    Multimap<String, Protein> proteins = HashMultimap.create();
    for (String f : files)
    {
      FastaProteinReader fasta = new CustomFastaProteinReader(new File(f));
      while (fasta.hasNext())
      {
        Protein p = fasta.next();
        proteins.put(p.toSymbolString(), p);
      }
      fasta.close();
    }
    return proteins;
  }

}
