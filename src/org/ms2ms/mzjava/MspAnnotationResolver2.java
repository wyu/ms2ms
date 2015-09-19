package org.ms2ms.mzjava;

/**
 * Copyright (c) 2010, SIB. All rights reserved.
 *
 * SIB (Swiss Institute of Bioinformatics) - http://www.isb-sib.ch Host -
 * http://mzjava.expasy.org
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution. Neither the name of the SIB/GENEBIO nor the names of
 * its contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL SIB/GENEBIO BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

  import com.google.common.base.Optional;
  import org.expasy.mzjava.core.mol.Mass;
  import org.expasy.mzjava.core.mol.NumericMass;
  import org.expasy.mzjava.core.ms.spectrum.IonType;
  import org.expasy.mzjava.proteomics.io.ms.spectrum.sptxt.AnnotationResolver;
  import org.expasy.mzjava.proteomics.mol.Peptide;
  import org.expasy.mzjava.proteomics.ms.spectrum.PepFragAnnotation;
  import org.expasy.mzjava.proteomics.ms.spectrum.PepLibPeakAnnotation;
  import org.expasy.mzjava.utils.RegexConstants;

  import java.util.ArrayList;
  import java.util.Collections;
  import java.util.List;
  import java.util.regex.Matcher;
  import java.util.regex.Pattern;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class MspAnnotationResolver2 implements AnnotationResolver {

  // 1: type, 2: res#, 3: nl, 4: charge,
  // 5: isotope, 6: isotope#, 7: mz diff
  private final String regex =
    "([abcxyzp?]|I[A-Z]{1,2})((?!-)" + RegexConstants.INTEGER
      + ")?((?=[-+])" + RegexConstants.INTEGER + ")?\\*?(i)?(\" + RegexConstants.INTEGER + \"|\\+)?(?:\\^("
      + RegexConstants.INTEGER + "))?/(" + RegexConstants.REAL + ")";

  private final Pattern annotPattern = Pattern.compile(regex);

  @Override
  public List<PepLibPeakAnnotation> resolveAnnotations(String value, Peptide peptide) {

    if(value.startsWith("?")) return Collections.emptyList();

    final String[] split = value.split("\\s+")[0].split(",");
    List<PepLibPeakAnnotation> annotationList = new ArrayList<PepLibPeakAnnotation>(split.length);

    for (final String annotation : split) {
      PepLibPeakAnnotation peakAnnotation;

      // do not parse unknown annots and internals
      if (annotation.charAt(0) != '?' && !"Int".equals(annotation.substring(0, 3)) && annotation.charAt(0) != 'I') {

        peakAnnotation = new PepLibPeakAnnotation(0, 0, 0, Optional.fromNullable(parseAnnotation(annotation, peptide)));

        annotationList.add(peakAnnotation);
      }
    }

    return annotationList;
  }

  /**
   * Parse annotation from the given string.
   *
   * @param annot the annotation to parse.
   * @return an annotation.
   */ //todo this should read the stats and return a LibraryPeakAnnotations
  public PepFragAnnotation parseAnnotation(final String annot, final Peptide peptide) {

    IonType ionType = null;
    int residueNumber = 0;
    double neutralLossMass = 0;
    int charge = 1;
    int c13Number = 0;

    final Matcher matcher = annotPattern.matcher(annot);

    if (matcher.find()) {
      /** 1. Get fragment type */
      final char typeChar = matcher.group(1).charAt(0);

      // Immonium
      if (typeChar == 'I') {
        final char residueChar = matcher.group(1).charAt(1);

        ionType = IonType.i;
        residueNumber = 1;
      } else {
        switch (typeChar) {
          case 'a':
            ionType = IonType.a;
            break;
          case 'b':
            ionType = IonType.b;
            break;
          case 'c':
            ionType = IonType.c;
            break;
          case 'x':
            ionType = IonType.x;
            break;
          case 'y':
            ionType = IonType.y;
            break;
          case 'z':
            ionType = IonType.z;
            break;
          case 'p':
            // parental
            ionType = IonType.p;
            residueNumber = peptide.size();
            break;
          case '?':
            // no annotation for this peak
            return null;
          default:
            throw new IllegalStateException(typeChar + " : bad ion type ");
        }

        /** 2. Get number of residues */
        if (matcher.group(2) != null) {
          residueNumber = Integer.parseInt(matcher.group(2));
        }
      }

      /** 3. Get neutral loss (da) */
      if (matcher.group(3) != null) {
        neutralLossMass = Double.parseDouble(matcher.group(3));
      }

      /** 4. C13 inside */
      if (matcher.group(4) != null) {

        c13Number = 1;
        /**
         * 6. number of C13 (with optional number of c13 (int) or '+' if
         * unknown)
         */

        if (matcher.group(5) != null) {
          if (matcher.group(5).length() != 0) {
            // if c13Number not known, use 1 for the moment
            // better to set it to the closest c13Number
            if (!"+".equals(matcher.group(5))) {
              c13Number = Integer.parseInt(matcher.group(6));
            }
          }
        }
      }

      /** 6. charge number */
      if (matcher.group(6) != null) {
        charge = Integer.parseInt(matcher.group(6));
      }

      /** 7. exp-theo peak mz */
    }
    if (ionType == null) {
      return null;
    }

    Peptide peptideFragment;
    switch (ionType.getFragmentType()) {

      case FORWARD:
      case REVERSE:
        peptideFragment = peptide.subSequence(ionType, residueNumber);
        break;
      case MONOMER:
        peptideFragment = peptide.subSequence(residueNumber-1,residueNumber);
        break;
      case INTACT:

        peptideFragment = new Peptide(peptide);
        break;
      case INTERNAL:
      default:
        throw new IllegalArgumentException("Cannot create annotation for ions of type " + ionType);
    }

    Mass neutralLoss = Mass.ZERO;
    if (neutralLossMass < 0) {

      neutralLoss = new NumericMass(neutralLossMass);
    }

    return new PepFragAnnotation.Builder(ionType, charge, peptideFragment).addC13(c13Number).setNeutralLoss(neutralLoss).build();
  }

}
