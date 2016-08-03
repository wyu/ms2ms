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

import java.util.*;

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
  public static String keepAAs(String seq)
  {
    String peptide=null;
    if (Strs.isSet(seq))
    {
      for (char c : seq.toCharArray())
        if (sAA.indexOf(c)>=0) peptide=Strs.extend(peptide, c+"","");
    }
    return peptide;
  }
  // setup a simple mapping of the AA symbol and their incremental masses
  public static Map<Character, Double> newAAsMass(String fixed)
  {
    Map<Character, Double> AAs = new HashMap<>();
    AAs.put('G', 57.02146d );AAs.put('A', 71.03711d );AAs.put('S', 87.03203d );AAs.put('P', 97.05276d );AAs.put('V', 99.06841d);
    AAs.put('T', 101.04768d);AAs.put('L', 113.08406d);AAs.put('N', 114.04293d);AAs.put('D', 115.02694d);AAs.put('Q', 128.05858d);
    AAs.put('E', 129.04259d);AAs.put('M', 131.04049d);AAs.put('H', 137.05891d);AAs.put('F', 147.06841d);AAs.put('R', 156.10111d);
    AAs.put('C', 160.03065d);AAs.put('Y', 163.06333d);AAs.put('W', 186.07931d);AAs.put('K', 128.09496d);
    AAs.put('^', 1.00783d  ); // N-terminal mod
    AAs.put('$', 17.00273d ); // C-terminal mod

    if ("tmt10".equalsIgnoreCase(fixed) || "tmt6".equalsIgnoreCase(fixed))
    {
      AAs.put('K', 357.25789d);
      AAs.put('^', 229.162932d); // TMT-10
    }
    return AAs;
  }
}
