package org.ms2ms.data.ms;

import org.jgrapht.ext.EdgeNameProvider;

/**
 * Created by yuw on 9/15/16.
 */
public class DenovoEdgeLabelProvider<E extends DenovoEdge> implements EdgeNameProvider<E>
{
 public DenovoEdgeLabelProvider() { }

  /**
   * Returns the String representation an edge.
   *
   * @param edge the edge to be named
   */
  @Override public String getEdgeName(E edge)
  {
    return edge.getLabel();
  }
}
