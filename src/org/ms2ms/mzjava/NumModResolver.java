package org.ms2ms.mzjava;

import com.google.common.base.Optional;
import org.expasy.mzjava.core.mol.NumericMass;
import org.expasy.mzjava.core.mol.Weighable;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.mol.modification.NumericModificationResolver;
import org.expasy.mzjava.proteomics.mol.modification.unimod.UnimodModificationResolver;
import org.expasy.mzjava.proteomics.ms.ident.ModificationMatch;
import org.ms2ms.utils.Strs;

import java.util.Collections;
import java.util.Comparator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 9/28/14
 * Time: 8:31 PM
 * To change this template use File | Settings | File Templates.
 */
public class NumModResolver extends NumericModificationResolver
{
  private Pattern mModPattern = Pattern.compile("\\{?([-+]?\\d*\\.?\\d+)\\}?");
  private UnimodModificationResolver mUnimod = new UnimodModificationResolver();

  public NumModResolver()          { super(); }
  public NumModResolver(Pattern p) { super(); mModPattern = p; }
  public NumModResolver(String  s) { super(); mModPattern = Pattern.compile(s); }

  @Override
  public Optional<Modification> resolve(String input)
  {
    // short circuit a mod check bug for MSGF+ generated mzid file
    if (Strs.equals(input, "unknown modification")) return Optional.of(new Modification(input, new NumericMass(0d)));

    Matcher matcher = mModPattern.matcher(input);
    if (matcher.matches()) {

      double val = Double.parseDouble(matcher.group(1));
      return Optional.of(new Modification(input, new NumericMass(val)));
    }
    else
    {
      // default to Unimod entries
      Optional<Modification> r = mUnimod.resolve(input);
      if (r.isPresent()) return r;
    }

    return Optional.absent();
  }
}
