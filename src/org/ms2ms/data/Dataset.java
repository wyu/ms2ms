package org.ms2ms.data;

/** A logical collection of data that will be processed and produce a HData object for future analysis
 *
 * Author: wyu
 * Date: 10/3/14
 */
public interface Dataset
{
  public String getName();
  public HData  getData();
}
