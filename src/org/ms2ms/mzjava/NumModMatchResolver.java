package org.ms2ms.mzjava;

import com.google.common.base.Optional;
import org.expasy.mzjava.core.mol.NumericMass;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.ms.ident.ModificationMatch;
import org.expasy.mzjava.proteomics.ms.ident.ModificationMatchResolver;

/**
 * Created by yuw on 10/19/2015.
 */
public class NumModMatchResolver implements ModificationMatchResolver
{
  @Override
  public Optional<Modification> resolve(ModificationMatch modMatch)
  {
    double mass = modMatch.getMassShift();

    if (mass!=0)
      return Optional.of(new Modification(modMatch.getResidue().getSymbol()+(mass<0?"-":"")+mass, new NumericMass(mass)));

    return Optional.absent();
  }
}
