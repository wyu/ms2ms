package org.ms2ms.nosql;

import com.thinkaurelius.titan.core.TitanFactory;
import com.thinkaurelius.titan.core.TitanGraph;
import org.apache.commons.configuration.BaseConfiguration;
import org.apache.commons.configuration.Configuration;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 6/1/14
 * Time: 12:05 AM
 * To change this template use File | Settings | File Templates.
 */
abstract public class Titans
{
  public static TitanGraph openHBaseGraph()
  {
    Configuration conf = new BaseConfiguration();
    conf.setProperty("storage.backend","hbase");
    conf.setProperty("storage.directory","/media/data/titan");

    return TitanFactory.open(conf);
  }
  public static TitanGraph openGraph(String cfg)
  {
    // graph.properties. possibility for multiple graphs
    return TitanFactory.open(cfg);
  }
  public static void addPeptides()
  {
//    GeneNode geneNode = new GeneNode(manager.createNode(GeneNode.NODE_TYPE));
//    geneNode.setPositions(genePositionsSt);
//    graph.addEdge(null, genomeElementNode.getNode(), geneNode.getNode(), GenomeElementGeneRel.NAME);

  }
}
