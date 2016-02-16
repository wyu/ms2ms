package org.ms2ms.mzjava;

import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import org.ms2ms.utils.Tools;

import javax.annotation.Nonnull;

/** Not working!!!
 * Created by yuw on 10/22/2015.
 */
@Deprecated
public class SpectrumIdentifierByRunScan extends SpectrumIdentifier implements Comparable<SpectrumIdentifierByRunScan>
{
  public SpectrumIdentifierByRunScan(String spectrum)
  {
    super(spectrum);
  }
  public SpectrumIdentifierByRunScan(@Nonnull SpectrumIdentifier id)
  {
    super(id.getSpectrum());

    // copy the contents
    if (id.getName(                ).isPresent()) setName(                id.getName().get());
    if (id.getPrecursorNeutralMass().isPresent()) setPrecursorNeutralMass(id.getPrecursorNeutralMass().get());
    if (id.getPrecursorIntensity(  ).isPresent()) setPrecursorIntensity(  id.getPrecursorIntensity().get());
    if (id.getPrecursorMz(         ).isPresent()) setPrecursorMz(         id.getPrecursorMz().get());
    if (id.getAssumedCharge(       ).isPresent()) setAssumedCharge(       id.getAssumedCharge().get());
    if (id.getIndex(               ).isPresent()) setIndex(               id.getIndex().get());
    if (id.getSpectrumFile(        ).isPresent()) setSpectrumFile(        id.getSpectrumFile().get());

    addScanNumber(   id.getScanNumbers());
    addRetentionTime(id.getRetentionTimes());
  }

  @Override
  public int hashCode()
  {
    return 31 * getScanNumbers().hashCode() + getSpectrumFile().hashCode();
  }

  @Override
  public int compareTo(SpectrumIdentifierByRunScan o)
  {
    int r = 0;
    if (getSpectrumFile().isPresent() && o.getSpectrumFile().isPresent())
      r = getSpectrumFile().get().compareTo(o.getSpectrumFile().get());
//    if (r==0 && Tools.isSet(getScanNumbers()) && Tools.isSet(o.getScanNumbers()))
//      r = Integer.compare(getScanNumbers().getFirst().getValue(), o.getScanNumbers().getFirst().getValue());

    return r;
  }
}
