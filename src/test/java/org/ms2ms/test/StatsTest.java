package org.ms2ms.test;

import com.google.common.collect.Lists;
import org.junit.Test;
import org.ms2ms.algo.MsStats;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Tools;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 7/17/14
 * Time: 6:43 AM
 * To change this template use File | Settings | File Templates.
 */
public class StatsTest extends TestAbstract
{
  @Test
  public void hypergeometric() throws Exception
  {

  }
  @Test
  public void binomial_approximation() throws Exception
  {
    System.out.println("prob(20, 0.3, 100) >>> 0.016462853241869437\n1-prob(40-1, 0.3, 100) >>> 0.020988576003924564\n" +
            "normal_estimate(20, 0.3, 100) >>> 0.014548164531920815\nnormal_estimate 1-prob(40-1/100/0.3 >>> 0.024767304545069813)\n" +
        "poisson(20, 0.3, 100) >>> 0.013411150012837811\npoisson 1-prob(40-1/100/0.3) >>> 0.046253037645840323\n");

//  http://stackoverflow.com/questions/1095650/how-can-i-efficiently-calculate-the-binomial-cumulative-distribution-function
//  prob(20, 0.3, 100) >>> 0.016462853241869437

    System.out.println("binomial exact (20/100/0.3): "+Stats.binomial_exact(20, 100, 0.3));
    System.out.println("binomial_ln    (20/100/0.3): "+Stats.binomial(20, 100, 0.3));

//  1-prob(40-1, 0.3, 100) >>> 0.020988576003924564
    System.out.println("binomial exact 1-prob(40-1/100/0.3): "+(1-Stats.binomial_exact(40-1, 100, 0.3)));
    System.out.println("binomial_ln    1-prob(40-1/100/0.3): "+(1-Stats.binomial(      40-1, 100, 0.3)));

//  normal_estimate(20, 0.3, 100) >>> 0.014548164531920815
    System.out.println("binomial_normal_apprx (20/100/0.3): "+Stats.binomial_normal_estimate(20, 100, 0.3));
    System.out.println("binomial_normal_apprx 1-prob(40-1/100/0.3): "+(1-Stats.binomial_normal_estimate(40-1, 100, 0.3)));

//  poisson(20, 0.3, 100) >>> 0.013411150012837811
    System.out.println("binomial_poisson_apprx(20/100/0.3): "+Stats.binomial_poisson_estimate(20, 100, 0.3));
    System.out.println("binomial_poisson_apprx 1-prob(40-1/100/0.3): "+(1-Stats.binomial_poisson_estimate(40-1, 100, 0.3)));

  }
  @Test
  public void aggregator() throws Exception
  {
    System.out.println(MsStats.aggregate(Lists.newArrayList(2.3, 5.6, 2.5), MsStats.Aggregator.COUNT));
    System.out.println(MsStats.aggregate(Lists.newArrayList(2.3, "5.6", 2.5), MsStats.Aggregator.COUNT));
    System.out.println(MsStats.aggregate(Lists.newArrayList(2.3, 5.6, 2.5), MsStats.Aggregator.MEAN));
    System.out.println(MsStats.aggregate(Lists.newArrayList(2.3, "5.6", 2.5), MsStats.Aggregator.MEAN));
    System.out.println(MsStats.aggregate(Lists.newArrayList(2.3, "5.6TT", 2.5), MsStats.Aggregator.MEAN));
  }
  @Test
  public void interpolating() throws Exception
  {
    double[] xs = { 12.024,14.841,16.985,17.272,17.508,17.512,18.168,18.34,18.751,19.848,20.436,20.585,21.246,21.342,21.683,22.41,23.082,23.342,23.453,23.497,24.049,24.755,24.88,25.42 },
             ys = { 1.758,0.66049,0.77216,0.32907,0.68285,0.86646,0.53094,0.37996,0.80291,1.2129,0.94221,0.65332,1.8598,0.56271,3.4555,3.0465,0.92357,2.3241,6.443,0.69866,0.47629,0.37924,0.68876,1.0595 },
             Xs = {10.933,14.432,16.669,17.067,17.303,17.117,17.885,18.128,18.461,19.196,19.931,20.215,20.345,21.065,20.596,21.073,22.536,22.388,21.045,23.061,23.766,24.531,24.641,24.81 };

    double[] Ys = MsStats.interpolate(xs, ys, 0.3, Xs);

    for (int i=0; i<xs.length; i++)
    {
      System.out.println(i + "\t" + Tools.d2s(xs[i], 3) + "\t" + Tools.d2s(ys[i], 3) + "\t"+ Tools.d2s(Xs[i], 3) + "\t" + Tools.d2s(Ys[i], 3));
    }
    System.out.println();
  }
}
