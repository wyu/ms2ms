package org.ms2ms.mzjava;

import com.google.common.base.Optional;
import org.expasy.mzjava.core.mol.NumericMass;
import org.expasy.mzjava.proteomics.io.ms.ident.ProteinPilotModificationResolver;
import org.expasy.mzjava.proteomics.mol.modification.Modification;

/**
 * Created by yuw on 11/13/2015.
 */
public class ProteinPilotModResolver extends ProteinPilotModificationResolver
{
  public ProteinPilotModResolver()
  {
    super();

    putOverrideUnimod("Chloro",             new Modification("Chloro",             new NumericMass( 33.961028d)));
    putOverrideUnimod("No Carbamidomethyl", new Modification("No Carbamidomethyl", new NumericMass(-57.021464d)));
    putOverrideUnimod("Methyl+Propionyl",   new Modification("Methyl+Propionyl",   new NumericMass( 70.041865d)));
    putOverrideUnimod("Phosphoglyceryl",    new Modification("Phosphoglyceryl",    new NumericMass(167.982375d)));
    putOverrideUnimod("Dichloro",           new Modification("Dichloro",           new NumericMass( 67.922055d)));
    putOverrideUnimod("Sec->Dha",           new Modification("Sec->Dha",           new NumericMass(-81.932171d)));
    putOverrideUnimod("No TMT6plex",        new Modification("No TMT6plex",        new NumericMass(-229.162932)));


//    putTranslate("Protein Terminal Acetyl", "Acetyl");

    // Transpeptidation
//    putOverrideUnimod("Arg-add", Modification.parseModification("Arg-add:C6H12N4O"));
//    putOverrideUnimod("Arg-loss", Modification.parseModification("Arg-loss:C-6H-12N-4O-1"));
//    putOverrideUnimod("No Methylthio", Modification.parseModification("No Methylthio:H-2C-1S-1"));
//    putOverrideUnimod("GlnGlnGlnThrGlyGly", Modification.parseModification("QQQTGG:H37C23N9O10")); //H(37) C(23) N(9) O(10)
//    putOverrideUnimod("Met->Hcy", new Modification("Met->Hcy", new NumericMass(-14.01565006)));

    putOverrideUnimod("Lys-add", Modification.parseModification("Lys-add:C6H12N2O"));
    putOverrideUnimod("GG",      Modification.parseModification("GG:H6C4N2O2"));
  }

  @Override
  public Optional<Modification> resolve(String input)
  {
    Optional<Modification> mod = super.resolve(input);

    // do something if no match is found
    if (!mod.isPresent())
    {

    }
    return mod;
  }

}
