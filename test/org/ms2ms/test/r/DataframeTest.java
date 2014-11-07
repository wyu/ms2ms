package org.ms2ms.test.r;

import org.expasy.mzjava.core.ms.PpmTolerance;
import org.expasy.mzjava.core.ms.Tolerance;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.ms2ms.r.Dataframe;
import org.ms2ms.runner.Aligner;
import org.ms2ms.test.TestAbstract;
import org.ms2ms.utils.Stats;

import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 7/17/14
 * Time: 6:42 AM
 * To change this template use File | Settings | File Templates.
 */
public class DataframeTest extends TestAbstract
{
  Dataframe obs1, obs2, animals;

  @Before
  public void setUp()
  {
    animals   = Dataframe.readtable(true, "size   type   name\nsmall  cat    lynx\nbig    cat    tiger\nsmall  dog    chihuahua\nbig    dog   greatdane");
    obs1      = Dataframe.readtable(true, "number size  type\n1      big   cat\n2      small dog\n3      small dog\n4      big   dog");
    obs2      = Dataframe.readtable(true, "number  size    type\n1       big     cat\n2       small   dog\n3       small   dog\n4       big     dog\n5       big     dog\n6       big     dog");
  }

  /**
   * animals *
   size   type   name
   small  cat    lynx
   big    cat    tiger
   small  dog    chihuahua
   big    dog   "great dane"

   * observations *
   number size  type
   1      big   cat
   2      small dog
   3      small dog
   4      big   dog

   * obs2 *
   number  size    type
   1       big     cat
   2       small   dog
   3       small   dog
   4       big     dog
   5       big     dog
   6       big     dog

   merge(observations, animals, c("size","type"))
   size   type  number    name
   big    cat   1         tiger
   big    dog   4         great dane
   small  dog   2         chihuahua
   small  dog   3         chihuahua

   > merge(obs2, animals, "size")
   size    number  type.x  type.y    name
   1      big     1       cat     cat       tiger
   2      big     1       cat     dog       great dane
   3      big     4       dog     cat       tiger
   4      big     4       dog     dog       great dane
   5      big     5       dog     cat       tiger
   6      big     5       dog     dog       great dane
   7      big     6       dog     cat       tiger
   8      big     6       dog     dog       great dane
   9      small   2       dog     cat       lynx
   10     small   2       dog     dog       chihuahua
   11     small   3       dog     cat       lynx
   12     small   3       dog     dog       chihuahua

   > merge(animals, obs2, "size")
   size   type.x  name         number  type.y
   1      big     cat   tiger         1       cat
   2      big     cat   tiger         4       dog
   3      big     cat   tiger         5       dog
   4      big     cat   tiger         6       dog
   5      big     dog   great dane    1       cat
   6      big     dog   great dane    4       dog
   7      big     dog   great dane    5       dog
   8      big     dog   great dane    6       dog
   9      small   cat   lynx          2       dog
   10     small   cat   lynx          3       dog
   11     small   dog   chihuahua     2       dog
   12     small   dog   chihuahua     3       dog
   *
   > merge(animals, obs2, by=c("size", "type"))
   size type       name number
   1   big  cat      tiger      1
   2   big  dog great dane      4
   3   big  dog great dane      5
   4   big  dog great dane      6
   5 small  dog  chihuahua      2
   6 small  dog  chihuahua      3
   *
   * @throws Exception
   */
  @Test
  public void merging() throws Exception
  {
    System.out.println(animals.display());
    System.out.println(obs1.display());
    System.out.println(obs2.display());

    Dataframe m1 = Dataframe.merge(animals, obs1, false);
    Dataframe m2 = Dataframe.merge(animals, obs2, false);
    Dataframe m3 = Dataframe.merge(animals, obs2, false, "size");

    System.out.println(m1.display());
    System.out.println(m2.display());
    System.out.println(m3.display());

    Assert.assertEquals(m1.rows().size(), 4);
    Assert.assertEquals(m2.rows().size(), 12);
    Assert.assertEquals(m3.rows().size(), 8);
  }

  /**
   * pivot( [dose sbj], visit_name ) produces the following table

   []               []    'visit_name'    'visit_name'
   'dose'      'sbj'            'D0'            'D22'
   'dosed'    '1003'    [         1]    [         1]
   'dosed'    '1015'    [         1]    [         1]
   'dosed'    '1025'    [         1]    [         1]
   *
   col is a categorical column whose factors will be used as the column header in the outgoing data frame
   val is a numberic column whose values will be the cell in the outgoing data frame
   func is the aggregate function if multiple values are found in a cell
   rows are the columns that will transferred to the outgoing data frame
   */
  @Test
  public void pivoting() throws Exception
  {
    Dataframe out = animals.pivot("type", "size", Stats.Aggregator.COUNT);
    System.out.println("\n" + out.display());
  }

  /** split(animals, animals$type)
   $cat
   size type  name
   1 small  cat  lynx
   2   big  cat tiger

   $dog
   size type       name
   3 small  dog  chihuahua
   4   big  dog great dane
   */
  @Test
  public void splitting() throws Exception
  {
    Map<Object, Dataframe> outs = animals.split("type");

    for (Object obj : outs.keySet())
      System.out.println(obj.toString() + "\n" + outs.get(obj).display());
  }
  /** ftable(animals)
   *
   name chihuahua great dane lynx tiger
   size  type
   big   cat               0          0    0     1
   dog               0          1    0     0
   small cat               0          0    1     0
   dog               1          0    0     0   *
   *
   */
  @Test
  public void ftabling() throws Exception
  {

  }

}
