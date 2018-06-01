package org.ms2ms.algo;

import libsvm.*;
import org.ms2ms.utils.Tools;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

// prediction of peptide retention time
//
public class PeptideLC
{
  public static void LsRP()
  {
//    Support Vector Regression. The ε support vector regression (ε-SVR) functionality in the libsvm package
//    was used in the RT predictor. A radial basis function kernel was employed and optimization of SVR parameters
//    is crucial to avoid overfitted or underfitted models. The best appropriate values of the parameters C, γ and ε
//    were figured out via 3-fold internal cross-validation, with C ∈ {2i | i∈ {−8,−7,…,8}}, γ ∈ {2i |i∈ {−8,−7,…,8}},
//    ε ∈ {10i |i∈ {−3,−2,−1}}.
  }
  // translate the peptide sequence into 500-cols vector
  public static svm_node[] Peptide2NN(String peptide)
  {
    if (peptide.length()>25) return null;

    svm_node[] nodes = newSVMnodes(500, 0d);
    // locate the mid point
    int mid = (int )Math.round(peptide.length()/2d), offset = 25-peptide.length();
    for (int k=0; k<peptide.length(); k++)
    {
      int pos = Peptides.sAA.indexOf(peptide.charAt(k));
      if (pos<0) return null;

      nodes[20*(k<mid?k:(k+offset))+pos].value = 1;
    }

    return nodes;
  }
  private static svm_node[] newSVMnodes(int size, double val)
  {
    svm_node[] nodes = new svm_node[size];
    for (int i=0; i<size; i++)
    {
      nodes[i] = new svm_node();
      nodes[i].index=i+1;
      nodes[i].value=val;
    }
    return nodes;
  }
  // taken from svm_toy.java in the distribution
  public static svm_parameter newSVRparam(double C, double gamma, double epsilon)
  {
    svm_parameter param = new svm_parameter();

    // default values
    param.svm_type = svm_parameter.EPSILON_SVR;
    param.kernel_type = svm_parameter.RBF;
    param.degree = 3;
//  param.C=1; param.gamma=1; param.eps=1e-3;
    param.C = C;
    param.gamma = gamma;
    param.eps = epsilon;
    param.coef0 = 0;
    param.nu = 0.5;
    param.cache_size = 40;
    param.p = 0.1;
    param.shrinking = 1;
    param.probability = 0;
    param.nr_weight = 0;
    param.weight_label = new int[0];
    param.weight = new double[0];

    return param;
  }

  // split the peptides into train and test according to the ratio
  public static void SVR(Map<String, Float> peptide_rt, float split)
  {
    // turn off the internal messages
    svm_print_interface print_func = svm_print_null;

    svm_parameter param = newSVRparam(Math.pow(2,6),Math.pow(2,0),Math.pow(10,-4));

    // build problem
    List<svm_node[]> x = new ArrayList<>(), x1 = new ArrayList<>();
    List<Float>      y = new ArrayList<>(), y1 = new ArrayList<>();

    Random RND = new Random(System.nanoTime());
    for(String peptide : peptide_rt.keySet())
    {
      svm_node[] nodes = Peptide2NN(peptide);
      if (nodes!=null)
      {
        if (RND.nextDouble()<split) { x.add(nodes);  y.add(peptide_rt.get(peptide)); }
        else                       { x1.add(nodes); y1.add(peptide_rt.get(peptide)); }
      }
    }
    if (!Tools.isSet(x) || !Tools.isSet(y) || x.size()!=y.size()) return;

    svm_problem prob = new svm_problem(), test = new svm_problem();
    prob.l = x.size();
    prob.y = new double[prob.l];
    prob.x = new svm_node[prob.l][500];

    for (int i=0; i<x.size(); i++)
    {
      prob.x[i] = x.get(i); prob.y[i] = y.get(i);
    }

    test.l = x1.size();
    test.y = new double[test.l];
    test.x = new svm_node[test.l][500];

    for (int i=0; i<x1.size(); i++)
    {
      test.x[i] = x1.get(i); test.y[i] = y1.get(i);
    }

    // build model & classify
    try
    {
      screenSVRparam(prob, param, 3);
//      svm_model model = svm.svm_train(prob, param);
//      xvalidateSVR(prob, param, 3);
      // check the performance!
//      System.out.println(model.toString());

//      predictSVR(model, prob);
//      predictSVR(model, test);
    }
    catch (Exception e)
    {
      e.printStackTrace();
    }

  }

  public static void predictSVR(svm_model model, svm_problem prob)
  {
    System.out.println("Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma="+svm.svm_get_svr_probability(model)+"\n");

    for (int i=0; i<prob.l; i++)
    {
      double v = svm.svm_predict(model,prob.x[i]);
      System.out.println(Tools.d2s(prob.y[i], 2)+"\t"+Tools.d2s(v, 2));
    }
  }
  private static void screenSVRparam(svm_problem prob, svm_parameter param, int nr_fold)
  {
    // default param first
////    double[] best = xvalidateSVR(prob, param, nr_fold);
//    double[] best = insituSVR(prob, param);
//    System.out.println("Screening @ C="+param.C+", gamma="+param.gamma+", epsilon="+param.eps);
//    System.out.print("Cross Validation Mean squared error = "+best[0]+"\n");
//    System.out.print("Cross Validation Squared correlation coefficient = "+ best[1]+"\n\n");

    // with C ∈ {2i | i∈ {−8,−7,…,8}}, γ ∈ {2i |i∈ {−8,−7,…,8}}, ε ∈ {10i |i∈ {−3,−2,−1}}.
    double[] best=null; StringBuffer buf = new StringBuffer();
    for (double c0=-8; c0<=24; c0++)
      for (double g0=-16; g0<=8; g0++)
        for (double e0=-4; e0<=-1; e0++)
        {
          param.C = Math.pow(2d, c0);
          param.gamma = Math.pow(2d, g0);
          param.eps = Math.pow(10d, e0);
          // start the cross-validation
          double[] r = xvalidateSVR(prob, param, nr_fold);
//          double[] r = insituSVR(prob, param);
          if (best==null || (best[1]<r[1]))
          {
//            buf.append("Screening @ C="+param.C+", gamma="+param.gamma+", epsilon="+param.eps+"\n");
//            buf.append("Cross Validation Mean squared error = "+r[0]+"\n");
//            buf.append("Cross Validation Squared correlation coefficient = "+ r[1]+"\n");
            buf.append("param\t"+c0+"\t"+g0+"\t"+e0+"\t"+r[0]+"\n"+r[1]+"\n");
            best=r;
          }
        }

     System.out.println("\n\n"+buf.toString());
  }
  private static double[] insituSVR(svm_problem prob, svm_parameter param)
  {
    double total_error = 0, sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;

    svm_model model = svm.svm_train(prob, param);

    for(int i=0;i<prob.l;i++)
    {
      double y = prob.y[i];
      double v = svm.svm_predict(model,prob.x[i]);
      total_error += (v-y)*(v-y);
      sumv += v;
      sumy += y;
      sumvv += v*v;
      sumyy += y*y;
      sumvy += v*y;
    }
    double r2 = ((prob.l*sumvy-sumv*sumy)*(prob.l*sumvy-sumv*sumy))/
        ((prob.l*sumvv-sumv*sumv)*(prob.l*sumyy-sumy*sumy));

    return new double[] {total_error/prob.l, r2};
  }

  private static double[] xvalidateSVR(svm_problem prob, svm_parameter param, int nr_fold)
  {
    int i;
    double total_error = 0;
    double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;
    double[] target = new double[prob.l];

    svm.svm_cross_validation(prob,param,nr_fold,target);
    for(i=0;i<prob.l;i++)
    {
      double y = prob.y[i];
      double v = target[i];
      total_error += (v-y)*(v-y);
      sumv += v;
      sumy += y;
      sumvv += v*v;
      sumyy += y*y;
      sumvy += v*y;
    }
    double r2 = ((prob.l*sumvy-sumv*sumy)*(prob.l*sumvy-sumv*sumy))/
        ((prob.l*sumvv-sumv*sumv)*(prob.l*sumyy-sumy*sumy));

    return new double[] {total_error/prob.l, r2};
  }
  private static svm_print_interface svm_print_null = new svm_print_interface()
  {
    public void print(String s) {}
  };
}
