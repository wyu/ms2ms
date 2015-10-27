package org.ms2ms.algo;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.expasy.mzjava.core.mol.SymbolSequence;
import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.expasy.mzjava.proteomics.mol.Protein;
import org.expasy.mzjava.proteomics.mol.modification.ModAttachment;
import org.expasy.mzjava.proteomics.ms.fragment.PeptideFragmentAnnotator;
import org.expasy.mzjava.proteomics.ms.fragment.PeptideFragmenter;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.util.Collection;
import java.util.EnumSet;
import java.util.List;

/**
 * Created by yuw on 9/14/2015.
 */
public class Peptides
{
  public static final String sAA = "ACDEFGHIKLMNPQRSTVWY";
  public static <P extends SymbolSequence> Multimap<String, P> toSequenceMap(Collection<P>... proteins)
  {
    if (!Tools.isSet(proteins)) return null;
    // create a map
    Multimap<String, P> seq_prots = HashMultimap.create();
    for (Collection<P> ps : proteins)
      for (P p : ps) seq_prots.put(p.toSymbolString(), p);

    return seq_prots;
  }
  public static PeptideFragmentAnnotator newHCDAnnotator()
  {
    // the fragmenter is responsible for generating a theoretical spectrum from the peptide
    PeptideFragmenter fragmenter =
      new PeptideFragmenter(EnumSet.of(IonType.b, IonType.y,IonType.a), PeakList.Precision.DOUBLE);

    // the annotator needs to delegate to the fragmenter the fragmentation process and to
    // the internal aligner the tolerance for aligning peaks
    return new PeptideFragmentAnnotator(fragmenter, new AbsoluteTolerance(0.1));
  }
}
