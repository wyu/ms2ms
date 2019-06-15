package org.ms2ms.test.ms;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.junit.Test;
import org.ms2ms.data.collect.TreeListMultimap;
import org.ms2ms.data.ms.FragmentEntry;
import org.ms2ms.data.ms.FragmentMatch;
import org.ms2ms.data.ms.FragmentPile;
import org.ms2ms.test.TestAbstract;
import toools.collections.Collections;

import java.util.Collection;
import java.util.Random;

public class PileTest extends TestAbstract
{
  FragmentPile pile = new FragmentPile(100);

  @Test
  public void FragmentPileTestRepeat() throws Exception
  {
    for (int k=0; k<5; k++)
    {
      System.out.println("Round "+k);
      FragmentPileTest();
    }

  }

  @Test
  public void FragmentPileTest() throws Exception
  {
    Random rnd = new Random(System.nanoTime());
    TreeListMultimap<Integer, FragmentMatch> data = TreeListMultimap.create();
    Multimap<Integer, Integer> pile_key = HashMultimap.create();

    pile.init();
    for (int i=0; i<28; i++)
    {
      int key = rnd.nextInt(), block = rnd.nextInt(10);
      for (int j=0; j<block; j++)
      {
        FragmentMatch pm = new FragmentMatch(new FragmentEntry(0f, key, null, j), 1);
        pile.add(pm);
        data.put(key, pm);
        pile_key.put(pile.getCurrPile(), key);
      }
    }
    pile.sort();

//    for (Integer key : data.keySet())
//      System.out.println(key+" --> "+data.get(key).size());
//
//    System.out.println();

    // key unique to the pile 2
    Collection<Integer> problems = Collections.difference(pile_key.get(1), pile_key.get(0));

    pile.start(); int track_size=1, rows=0;
    while (track_size>0)
    {
      track_size=pile.nextChunk(); rows+=track_size;
      System.out.println(pile.getCurrentKey()+" --> "+track_size+"/"+rows);

      if (data.containsKey(pile.getCurrentKey()) &&
          data.get(pile.getCurrentKey()).size()==track_size)
      {
        data.remove(pile.getCurrentKey());
//        System.out.println(pile.getCurrentKey()+" --> "+track_size);
//        FpmEntry ff = new FpmEntry(pile.at(track_size-1).getEntry(), pile.getTrack(), track_size).inspect4screen(2d);
      }
    }
    System.out.println("data pile not recoverred: "+data.size());
  }
}
