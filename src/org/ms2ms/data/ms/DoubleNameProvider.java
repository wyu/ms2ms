package org.ms2ms.data.ms;

import org.jgrapht.event.GraphEdgeChangeEvent;
import org.jgrapht.event.GraphListener;
import org.jgrapht.ext.VertexNameProvider;
import org.ms2ms.utils.Tools;

/**
 * Created by yuw on 9/16/16.
 */
public class DoubleNameProvider<V> implements VertexNameProvider<V>
{
  private int mDecimal=2;

  public DoubleNameProvider(int decimal)
  {
    mDecimal=decimal;
  }

  @Override
  public String getName(V vertex) { return getVertexName(vertex); }
  /**
   * Returns the String representation of the unique integer representing a
   * vertex.
   *
   * @param vertex the vertex to be named
   *
   * @return the name of
   *
   * @see GraphListener#edgeAdded(GraphEdgeChangeEvent)
   */
  @Override public String getVertexName(V vertex)
  {
    return Tools.d2s(((Double )vertex), mDecimal);
  }
}