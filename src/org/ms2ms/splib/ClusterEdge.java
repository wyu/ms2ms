package org.ms2ms.splib;

import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.ms2ms.data.Binary;
import org.ms2ms.data.ms.Ms2Cluster;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

public class ClusterEdge extends DefaultWeightedEdge implements Binary
{
  public enum EdgeType { REF, SAME_MH, Z, C13, MOD, LM, EMBED, RAGGED, NONE };

//  private Ms2Cluster mSource, mTarget;
  private int mSharedScans=0, tSourceHash, tTargetHash;
  private float mForwardDP, mReverseDP, mDeltaMH;
  private String mLabel;
  private EdgeType mType = EdgeType.NONE;
  private static final long serialVersionUID = 3258409452177034555L;

  public ClusterEdge() { super(); }
  public String getLabel() { return mLabel; }
  public ClusterEdge setLabel(String s) { mLabel=s; return this; }

  public ClusterEdge setForwardDP(float s) { mForwardDP=s; return this; }
  public ClusterEdge setReverseDP(float s) { mReverseDP=s; return this; }
  public ClusterEdge setDeltaMH(  float s) { mDeltaMH=s; return this; }
  public ClusterEdge setSharedScans(int s) { mSharedScans=s; return this; }
  public ClusterEdge setType(EdgeType s) { mType=s; return this; }

  public EdgeType getType() { return mType; }
  public float getForwardDP() { return mForwardDP; }
  public float getReverseDP() { return mReverseDP; }
  public float getDeltaMH()   { return mDeltaMH; }

  @Override public String toString()
  {
    return getSource().toString() + "<"+(mType.equals(EdgeType.NONE)?"":mType.toString())+">" + getTarget().toString() + (Strs.isSet(mLabel)?"("+mLabel+")":"");
  }
  @Override
  public void write(DataOutput ds) throws IOException
  {
    IOs.write(ds, getSource()!=null?getSource().hashCode():null);
    IOs.write(ds, getTarget()!=null?getTarget().hashCode():null);
    IOs.write(ds, mSharedScans);
    IOs.write(ds, mLabel);
    IOs.write(ds, mType.name());
    IOs.write(ds, mForwardDP);
    IOs.write(ds, mReverseDP);
    IOs.write(ds, mDeltaMH);
  }

  @Override
  public void read(DataInput ds) throws IOException
  {
    tSourceHash = IOs.read(ds, 0);
    tTargetHash = IOs.read(ds, 0);
    mSharedScans= IOs.read(ds, 0);
    mLabel      = IOs.read(ds, mLabel);
    try
    {
      mType = EdgeType.valueOf(IOs.read(ds, "NONE"));
    }
    catch (IllegalArgumentException e)
    {
      mType = EdgeType.NONE;
    }
    mForwardDP = IOs.read(ds, mForwardDP);
    mReverseDP = IOs.read(ds, mReverseDP);
    mDeltaMH   = IOs.read(ds, mDeltaMH);
  }

  public static void write(DataOutput ds, Graph<Ms2Cluster, ClusterEdge> graph) throws IOException
  {
    if (graph==null || !Tools.isSet(graph.vertexSet())) { IOs.write(ds,0); return; }

    // still need a place holder
    IOs.write(ds,graph.vertexSet().size());

    Map<Integer, Ms2Cluster> hash_node = new HashMap<>();
    for (Ms2Cluster cls : graph.vertexSet()) hash_node.put(cls.hashCode(), cls);

    IOs.writeIntMap(ds, hash_node);
    IOs.write(ds, graph.edgeSet());
  }
  public static Graph<Ms2Cluster, ClusterEdge> read(DataInput ds, Graph<Ms2Cluster, ClusterEdge> graph) throws IOException
  {
    int n = IOs.read(ds, 0);
    if (n<=0) return graph;

    Map<Integer, Ms2Cluster> nodes = IOs.readIntMap(ds, new HashMap<>(), Ms2Cluster.class);
    Collection<ClusterEdge>  edges = IOs.read(ds, new ArrayList<>(),    ClusterEdge.class);

    // let's reconstruct the graph
    for (Ms2Cluster v : nodes.values()) graph.addVertex(v);
    for (ClusterEdge e : edges)
    {
      graph.addEdge(nodes.get(e.tSourceHash), nodes.get(e.tTargetHash), e);
    }
    return graph;
  }
}
