package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.peaklist.PeakAnnotation;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.expasy.mzjava.proteomics.ms.spectrum.LibrarySpectrum;
import org.ms2ms.Disposable;
import org.ms2ms.utils.Tools;

import java.io.Serializable;

/** Copyright 2014-2015 ms2ms.org
 * <p/>
 * Description:  Dedicated class to serve the binary I/O
 * <p/>
 * Author: wyu
 * Date:   11/13/14
 */
public class MsSpectrum  implements Serializable, Disposable
{
  private static final long serialVersionUID = 8472752523296641667L;

  protected int precursorZ;        // the precursor charge
  protected double[] mzList;            // m/z of the peaks
  protected short[] intensityList;     // relative intensities of the peaks
  protected float precursorAi;
  protected double precursorMz;         // the precursor charge
  protected int size;              // length of the peak list
  // the upper bound of the peak intensity and variance in m/z and intensity
  protected float maxIntensity;

  public <A extends PeakAnnotation> MsSpectrum(PeakList<A> src)
  {
    super();
    // TODO copy the data from the source
    precursorMz =        src.getPrecursor().getMz();
    precursorAi =(float )src.getPrecursor().getIntensity();
    precursorZ  =        src.getPrecursor().getCharge();

    mzList        = new double[src.size()];
    intensityList = new short[ src.size()];
    for (int i=0; i<src.size(); i++)
    {
      mzList[       i]=src.getMz(i);
      intensityList[i]=(short )(src.getIntensity(i)*Short.MAX_VALUE/maxIntensity+Short.MIN_VALUE);
    }
    maxIntensity=(float )src.getBasePeakIntensity();
    size        =        src.size();
  }

  /** compression and uncompression routines  **/
  protected double getMz(int i)
  {
    rangeCheck(i); return mzList[i];
  }

  protected double getIntensity(int i)
  {
    rangeCheck(i);
    return maxIntensity * (intensityList[i]-Short.MIN_VALUE)/(2d*Short.MAX_VALUE);
  }

  private void rangeCheck(int index)
  {
    if (index >= size) throw new IndexOutOfBoundsException("Index: " + index + ", Size: " + size);
  }

  @Override
  public void dispose()
  {
    intensityList=null; mzList=null;
  }
  public static MsSpectrum adopt(MsnSpectrum s) { return new MsSpectrum(s); }
}
