package org.ms2ms.test;

import com.thinkaurelius.titan.core.TitanGraph;
import com.tinkerpop.blueprints.Direction;
import com.tinkerpop.blueprints.Edge;
import com.tinkerpop.blueprints.Vertex;
import org.junit.Test;
import org.ms2ms.nosql.Titans;


/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 5/31/14
 * Time: 11:13 PM
 * To change this template use File | Settings | File Templates.
 */
public class TitanTest01 extends TestAbstract
{
  @Test
  public void BasicOps()
  {
    TitanGraph graph = Titans.openHBaseGraph();

    Vertex a = graph.addVertex(null);
    Vertex b = graph.addVertex(null);
    a.setProperty("name", "marko");
    b.setProperty("name", "peter");
    Edge e = graph.addEdge(null, a, b, "knows");
    System.out.println(e.getVertex(Direction.OUT).getProperty("name") + "--" + e.getLabel() + "-->" + e.getVertex(Direction.IN).getProperty("name"));
  }
}
