package com.amgen.seattle.bioinfo.algorithm;

/** Adapted from the Matlab scripts by
 *
 * Author:
 * Jianwen Luo <luojw@ieee.org> 2005-06-12
 * Department of Biomedical Engineering
 * Tsinghua University, Beijing 100084, P. R. China
 *
 * References:
 * [1]	Luo J W, Ying K, Bai J.
 *     Savitzky-Golay Smoothing and Differentiation Filter for Even Number Data,
 *     Signal Process. 2005, 85(7): 1429-1434
 * [2]	Luo J W, Bai J.
 *     Savitzky-Golay Smoothing and Differentiation Filter of Even Length: A Gram Polynomial Approach,
 *     Spectroscopy. 2005, in press
 *    (http://www.spectroscopymag.com/spectroscopy/article/articleDetail.jsp?id=197394)
 *
 * Author: Wen Yu
 * Date: Feb 8, 2007
 * Time: 7:17:21 AM
 * <p/>
 * Copyright &copy; 2004-2005 Immunex Corporation, an affiliate of Amgen Inc.
 */
public class sgsdf
{
  /**
   % Savitzky-Golay Smoothing and Differentiation Filter of Even Length
   % (Closed-Form Solution)
   %
   % n:      polynomial degree
   % s:      derivative(differentiation) order (0=smoothing)
   %
   % Author:
   % Jianwen Luo <luojw@ieee.org> 2005-06-12
   % Department of Biomedical Engineering
   % Tsinghua University, Beijing 100084, P. R. China
   %
   % References:
   % [1]	Luo J W, Ying K, Bai J.
   %     Savitzky-Golay Smoothing and Differentiation Filter for Even Number Data,
   %     Signal Process. 2005, 85(7): 1429-1434
   % [2]	Luo J W, Bai J.
   %     Savitzky-Golay Smoothing and Differentiation Filter of Even Length: A Gram Polynomial Approach,
   %     Spectroscopy. 2005, in press
   *
   * @param n
   * @param s
   */
  public static double[] sgsdf_even_closed_form(int n, int s)
  {
    /*
    i=sym('i');% -m+1<=i<=m
    m=sym('m');
    t=sym(0); %evaluated at the imaginary central point
    subs(factor(hinstm(i,n,s,t,m)),'i','i-1/2')
    */
    return new double[]{};
  }

  /** Savitzky-Golay Smoothing and Differentiation Filter of Even Length
   *
   * @param n =  polynomial degree
   * @param s =  derivative(differentiation) order (0=smoothing)
   * @param t =  evaluation point (commonly,t=0,i.e.,smoothing or differentiation
   *             at the imaginary central point of the 2*m points)
   * @param m = (2*m) => data point number,i.e., filter length
   * @return hnstm:  convoluction coefficients
   */
  public static double[] sgsdf_even_gram_poly(double n, double s, double t, double m)
  {
    double[] hnstm= new double[2*(int )m];
    for (int i=0; i < 2*m; i++)
      hnstm[i] = hinstm(i,n,s,t,m);

    return hnstm;
  }

  /** Equs. (8) and (9) in Ref.[2]
   *
   * @param i
   * @param s =  derivative(differentiation) order (0=smoothing)
   * @param t =  evaluation point (commonly,t=0,i.e.,smoothing or differentiation
   *             at the imaginary central point of the 2*m points)
   * @param m = (2*m) => data point number,i.e., filter length
   */
  public static double hinstm(double i, double n, double s, double t, double m)
  {
    double value = 0;
    for (int k=0; k < n; k++)
    {
      value += ((2*k+1) * gff(2*m-1,k) / gff(2*m+k,k+1)) * Pkmi(k,m,i) * Pkmsi(k,m,s,t);
    }
    return value;
  }
  /** Gram polynomials (Equ.(14) in Ref.[2])
   *
   * @param k
   * @param m = (2*m) => data point number,i.e., filter length
   * @param i
   */
  public static double Pkmi(double k, double m, double i)
  {
    double value = 0;
    if (k == -1) value = 0; else if (k == 0) value = 1;
    else
    {
      //value=2*(2*k-1)/k/(2*m-k)*i*Pkmi(k-1,m,i)-(k-1)*(2*m+k-1)/k/(2*m-k)*Pkmi(k-2,m,i);
      value= (2*(2*k-1)/(k*(2*m-k)))*i*Pkmi(k-1,m,i)-((k-1)*(2*m+k-1)/(k*(2*m-k)))*Pkmi(k-2,m,i);
    }
    return value;
  }
  /** s-th derivative of Gram polynomials (Equ. (15) in Ref.[2])
   *
   * @param k
   * @param m = (2*m) => data point number,i.e., filter length
   * @param s =  derivative(differentiation) order (0=smoothing)
   * @param i
   */
  public static double Pkmsi(double k, double m, double s, double i)
  {
    double value = 0;
    if (s==0) value = Pkmi(k,m,i); else if (k==-1 || k==0) value = 0;
    else
    {
      //value = 2*(2*k-1)/k/(2*m-k)*(i*Pkmsi(k-1,m,s,i)+s*Pkmsi(k-1,m,s-1,i))-(k-1)*(2*m+k-1)/k/(2*m-k)*Pkmsi(k-2,m,s,i);
      value = (2*(2*k-1)/(k*(2*m-k)))*(i*Pkmsi(k-1,m,s,i)+s*Pkmsi(k-1,m,s-1,i))-((k-1)*(2*m+k-1)/(k*(2*m-k)))*Pkmsi(k-2,m,s,i);
    }
    return value;
  }
  /** generalized factorial function (Equ.(10) in Ref.[2])
   *
   *    (b)   { 1 (b == 0)
   * (a)    = |
   *          { a(a-1)...(a-b+1) (b >= 1)
   */
  public static double gff(double a, double b)
  {
    double value = a;

    if (b == 0) value = 1;
    else
    {
      for (int B = 1; B < b; B++) value *= B;
    }
    return value;
  }
}
