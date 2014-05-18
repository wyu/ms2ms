package org.ms2ms.nosql;

import java.io.Serializable;

/** Parameters about the HBase table including name, row counts, primary key
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 5/14/14
 * Time: 8:36 PM
 * To change this template use File | Settings | File Templates.
 */
public class HBaseTable implements Serializable
{
  private String name;
  long idend, entries;

}
