package org.ms2ms.mzjava;

import org.expasy.mzjava.proteomics.io.ms.spectrum.sptxt.LibrarySpectrumBuilder;
import org.expasy.mzjava.proteomics.io.ms.spectrum.sptxt.SpectraStCommentParser;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 9/28/14
 * Time: 9:17 PM
 * To change this template use File | Settings | File Templates.
 */
public class SpCommentParser extends SpectraStCommentParser {

  @Override
  protected boolean parseUnknownCommentTag(String tag, String value, LibrarySpectrumBuilder builder)
  {
    // TODO need to customize the tags
//    if("Mz_exact".equals(tag)) parseMzExact(value, builder);
    return false;
  }

  private void parseMzExact(String value, LibrarySpectrumBuilder builder) {

    builder.setPrecursorMz(Double.parseDouble(value.trim()));
  }
}
