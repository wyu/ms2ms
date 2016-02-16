package org.ms2ms.test.ms;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import org.expasy.mzjava.core.mol.Composition;
import org.expasy.mzjava.core.mol.Mass;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.expasy.mzjava.proteomics.mol.AminoAcid;
import org.expasy.mzjava.proteomics.mol.ModifiedPeptideFactory;
import org.expasy.mzjava.proteomics.mol.Peptide;
import org.expasy.mzjava.proteomics.mol.Protein;
import org.expasy.mzjava.proteomics.mol.digest.Protease;
import org.expasy.mzjava.proteomics.mol.digest.ProteinDigester;
import org.expasy.mzjava.proteomics.mol.modification.ModAttachment;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.ms.fragment.BackbonePeakGenerator;
import org.expasy.mzjava.proteomics.ms.fragment.PeptideFragmenter;
import org.expasy.mzjava.proteomics.ms.fragment.PeptideNeutralLossPeakGenerator;
import org.expasy.mzjava.proteomics.ms.fragment.PeptidePeakGenerator;
import org.expasy.mzjava.proteomics.ms.spectrum.PepFragAnnotation;
import org.expasy.mzjava.proteomics.ms.spectrum.PeptideSpectrum;
import org.junit.Before;
import org.junit.Test;
import org.ms2ms.io.MsIO;
import org.ms2ms.math.Histogram;
import org.ms2ms.test.TestAbstract;
import org.ms2ms.utils.Tools;

import java.util.*;

/**
 * Created by yuw on 10/31/2015.
 */
public class MassDefectTest extends TestAbstract
{
  ProteinDigester     digester = null;
  PeptideFragmenter fragmenter = null;

  @Before
  public void setUp() throws Exception
  {
    digester = new ProteinDigester.Builder(Protease.TRYPSIN).build();

    Set<IonType> ionTypes = EnumSet.of(IonType.b,IonType.y,IonType.a, IonType.i, IonType.p);

    Mass    waterLoss = Composition.parseComposition("H-2O-1");
    Mass ammoniumLoss = Composition.parseComposition("H-3N-1");

    List<PeptidePeakGenerator<PepFragAnnotation>> peakGenerators = Lists.newArrayList();

    peakGenerators.add(new BackbonePeakGenerator(ionTypes, 1));
    peakGenerators.add(new PeptideNeutralLossPeakGenerator(waterLoss,
      EnumSet.of(AminoAcid.S, AminoAcid.T, AminoAcid.D, AminoAcid.E), ionTypes, 1));
    peakGenerators.add(new PeptideNeutralLossPeakGenerator(ammoniumLoss,
      EnumSet.of(AminoAcid.Q, AminoAcid.K, AminoAcid.R, AminoAcid.N), ionTypes, 1));

    fragmenter = new PeptideFragmenter(peakGenerators, PeakList.Precision.DOUBLE);
  }
  @Test
  public void MassDefectHCD() throws Exception
  {
    double step = 1d;
    Multimap<String, Protein> sequences = MsIO.readFASTA("C:/local/sequences/UPS1_UPS2_Ecoli_469008_RefSeq.fasta");

    // this factory generate set of ProteinDigestProduct given variable modifications
    ModifiedPeptideFactory factory = new ModifiedPeptideFactory.Builder().build();

    // we want a series of histograms by the fragment masses
    SortedMap<Long, Histogram> mass_defect = new TreeMap<>();
    // for each protein
    int counts=0;
    for (Protein protein : sequences.values())
    {
      List<Peptide> digests = digester.digest(protein);
      for (Peptide digested : digests)
      {
        if (digested.getMolecularMass()>3000) continue;

        PeptideSpectrum spectrum = fragmenter.fragment(digested, 1);
        for (int i = 0; i < spectrum.size(); i++)
        {
          long position = Math.round(spectrum.getMz(i)/step);
          Histogram hist = mass_defect.get(position);
          if (hist == null)
          {
            hist = new Histogram(50);
            mass_defect.put(position, hist);
          }
          // deposit the mass
          hist.add(spectrum.getMz(i) - Math.round(spectrum.getMz(i)));
        }
        spectrum.clear(); spectrum=null;
      }
      digests.clear(); digests=null;
      if (++counts%100 ==0) System.out.print(".");
      if (  counts%5000==0)
      {
        System.gc();
        System.out.println("$" + counts);
      }
    }
    System.out.println();

    // inspect the distribution of the mass defects
    for (Long position : mass_defect.keySet())
    {
      Histogram hist = mass_defect.get(position).generate();
      System.out.println(position*step +"\t"+ Tools.d2s(hist.getMean(), 4) +"\t"+
        Tools.d2s(1E6*hist.getMean() /(position*step), 4) + "\t" + Tools.d2s(hist.getStdev(), 4) + "\t" +
        Tools.d2s(1E9*hist.getStdev()/(position*step), 4) + "\t" + hist.getData().size());
    }
  }

}
