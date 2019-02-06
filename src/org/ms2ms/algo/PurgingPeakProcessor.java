package org.ms2ms.algo;

import org.expasy.mzjava.core.ms.peaklist.AbstractPeakProcessor;
import org.expasy.mzjava.core.ms.peaklist.PeakAnnotation;

import java.util.List;

/**
 * Created by yuw on 10/6/2015.
 */
public class PurgingPeakProcessor<A extends PeakAnnotation> extends AbstractPeakProcessor<A, A>
{
  @Override
  public void start(int size) { sink.start(size); }

  @Override
  public void processPeak(double mz, double intensity, List<A> annotations)
  {
    if (intensity>0) sink.processPeak(mz, intensity, annotations);
  }

  @Override
  public void end() { sink.end(); }
}
