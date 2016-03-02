package org.ms2ms.data.ms;

import com.compomics.util.experiment.identification.matches.PeptideMatch;
import com.google.common.collect.Multimap;

import java.util.Collection;

/**
 * Created by yuw on 2/24/16.
 */
public class ProteinID
{

  ProteinID mParent=null;
  Collection<ProteinID> mChildren=null;

  Multimap<String, PeptideMatch> mSeqMatch=null;
}
