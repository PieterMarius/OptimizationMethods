using System;
using MathNet.Numerics;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.RootFinding;
using ConsoleApplication1.Optimization;
using ConsoleApplication1.Optimization.SteepestDescentMethod;
using ConsoleApplication1.Optimization.NonLinearConjugateGradient;
using ConsoleApplication1.Optimization.LinearSystem;
using ConsoleApplication1.Optimization.SequentialQuadraticProgramming;

namespace ConsoleApplication1
{
    class Program
    {
        static void Main(string[] args)
        {
            double oneOver2Pi = 1.0 / (1.0 * Math.Sqrt(2 * Math.PI));
            
            Console.WriteLine(NormalCDFInverse(0.99));
            Console.WriteLine(inv_cdf(0.99));
            Console.WriteLine(InverseCDF.QNorm(0.99,0,1,true,false));
            Console.WriteLine(NormalCDFInverse(0.5));
            Console.WriteLine(inv_cdf(0.5));
            Console.WriteLine(InverseCDF.QNorm(0.5, 0, 1, true, false));
            Console.WriteLine(NormalCDFInverse(0.31));
            Console.WriteLine(inv_cdf(0.31));
            Console.WriteLine(InverseCDF.QNorm(0.31, 0, 1, true, false));

            //x4−8x2 + 5
            Func<Vector, double> testFunc1 = (x) =>
            {
                return Math.Pow(x.Vars[0], 4) - 8 * Math.Pow(x.Vars[0], 2) + 5;
            };

            Func<Vector, double>[] dtestFunc1 = new Func<Vector, double>[1];

            dtestFunc1[0] = (x) =>
            {
                return 4 * Math.Pow(x.Vars[0], 3) - 16 * x.Vars[0];
            };

            Func<Vector, double> testConstr = (x) =>
            {
                return 5 - Math.Exp(x.Vars[0]) + 2.0 * Math.Pow(x.Vars[0] - 1, 2);
            };

            Func<double[], double> testv = (x) =>
            {
                double a = x[0];
                return 0.0;
            };

            //Func<Variables, double> testFunc = (x) => {
            //    return Math.Pow(x.Vars[0], 4) - 3 * Math.Pow(x.Vars[0], 3) + 2;
            //};

            //Func<Variables, double> dTestfunc = (x) =>
            //{
            //    return 4 * Math.Pow(x.Vars[0], 3) - 9 * Math.Pow(x.Vars[0], 2);
            //};

            //x4−3x3 + 2

            //TestFunc();

            Func<Vector, double> bananaFunc = (x) =>
            {
                return Math.Pow(1 - x.Vars[0], 2) + 100 * Math.Pow(x.Vars[1] - x.Vars[0] * x.Vars[0], 2);
            };

            Func<Vector, double> powell = (x) =>
            {
                return Math.Pow(x.Vars[0] + 10*x.Vars[1], 2) + 
                       5 * Math.Pow(x.Vars[2] - x.Vars[3], 2) +
                       Math.Pow(x.Vars[1] + 2 * x.Vars[2], 4) +
                       10 * Math.Pow(x.Vars[0] - x.Vars[3], 4);
            };

            Vector[] ttt = new Vector[3];
            ttt[0] = new Vector(new double[3] { 1, 2 ,3 });
            ttt[1] = new Vector(new double[3] { 4, 5, 6 });
            ttt[2] = new Vector(new double[3] { 7, 8, 9 });

            var tr = Vector.Transpose(ttt);

            //TestsSQP.Test0();
            //TestsSQP.Test1();
            //TestsSQP.Test2();
            //TestsSQP.Test3();
            //TestsSQP.Test4();
            //TestsSQP.Test5();
            //TestsSQP.Test6();
            //TestsSQP.Test7();
            //TestsSQP.Test8();
            //TestsSQP.Test9();
            //TestsSQP.Test10();
            //TestsSQP.Test11();
            //TestsSQP.Test12();
            //TestsSQP.Test13();
            //TestsSQP.Test14();
            //TestsSQP.Test15();
            //TestsSQP.Test16();
            TestsSQP.Test17();
            //TestsSQP.Test18();
            TestsSQP.Test19();
            //TestCGMethod();

            BFGS bfsg = new BFGS();

            var  result0 = bfsg.Solve(powell, new double[] { 3.0, -1.0, 0.0, 1.0 }, 10000);

            Console.WriteLine("Min " + powell(result0));

            NLCG gradient = new NLCG();

            Vector result = gradient.Solve(powell,   new double[] { 3.0, -1.0, 0.0, 1.0 }, 4000);

            Console.WriteLine("Min " + powell(result));

            SteepestDescent gradientSteep = new SteepestDescent();

            Vector result1 = gradientSteep.Solve(powell, new double[] { 3.0, -1.0, 0.0, 1.0 }, 10000);

            Console.WriteLine("Min " + powell(result1));

            //Variables result = SteepestDescent(testFunc, dTestFunc, 2, 20);

            for (int i = 0; i < result.Vars.Length; i++)
                Console.WriteLine("result " + result.Vars[i]);

            //Console.WriteLine("ver " + testFunc(result));

            Console.ReadLine();
        }


        
        static void TestFunc()
        {
            Vector a = new Vector(new double[] { 1, 2 });
            Vector b = new Vector(new double[] { 3, 4 });

            Vector c = new Vector(new double[] { 1, 3 });
            Vector d = new Vector(new double[] { 4, 7 });

            var res = Vector.Mult(a, b);

            var res1 = Vector.SubtractFromIdentity(res);

            var res2 = Vector.Div(res, 2);

            Vector[] aa = new Vector[] { a, b };
            Vector[] bb = new Vector[] { c, d };

            var res3 = Vector.Mult(res, bb);
            var res4 = Vector.Sum(res, bb);

            var res5 = Vector.Mult(res, c);

            
        }

        

        static void TestCGMethod()
        {
            MINRES minres = new MINRES();
            CGMethod cg = new CGMethod();

            Vector[] A2 = new Vector[3];
            A2[0] = new Vector(new double[] { 2, 1, 3 });
            A2[1] = new Vector(new double[] { 2, 6, 8 });
            A2[2] = new Vector(new double[] { 6, 8, 18 });

            Vector b2 = new Vector(new double[] { 1, 3, 5 });

            var minresSol = minres.Solve(A2, b2, new Vector(new double[3]), 50);
            var sol2 = cg.Solve(A2, b2, new Vector(new double[3]), 50);

            Vector[] A = new Vector[3];
            A[0] = new Vector(new double[] { 2, 0, 1 });
            A[1] = new Vector(new double[] { 1, 6, 0 });
            A[2] = new Vector(new double[] { 3, 2, 3 });

            Vector x = new Vector(new double[] { 2, 5, 7 });
            
            var sol = cg.Solve(A, x, new Vector(new double[3]), 50);
            var solm = minres.Solve(A, x, new Vector(new double[3]), 50);

            Vector[] A1 = new Vector[3];
            A1[0] = new Vector(new double[] { 3, 1, -6 });
            A1[1] = new Vector(new double[] { 2, 1, -5 });
            A1[2] = new Vector(new double[] { 6, -3, 3 });

            Vector b1 = new Vector(new double[] { -10, -8, 0 });

            var sol1 = cg.Solve(A1, b1, new Vector(new double[3]), 200);
            var solm1 = minres.Solve(A1, b1, new Vector(new double[3]), 200);

            Vector diff = b1 - Vector.Mult(A1, solm1);
            Vector diff1 = b1 - Vector.Mult(A1, sol1);
        }

        static double RationalApproximation(double t)
        {
            // Abramowitz and Stegun formula 26.2.23.
            // The absolute value of the error should be less than 4.5 e-4.
            double[] c = { 2.515517, 0.802853, 0.010328 };
            double[] d = { 1.432788, 0.189269, 0.001308 };
            return t - ((c[2] * t + c[1]) * t + c[0]) /
                        (((d[2] * t + d[1]) * t + d[0]) * t + 1.0);
        }

        static double NormalCDFInverse(double p)
        {
            if (p <= 0.0 || p >= 1.0)
            {

            }

            // See article above for explanation of this section.
            if (p < 0.5)
            {
                // F^-1(p) = - G^-1(p)
                return -RationalApproximation(Math.Sqrt(-2.0 * Math.Log(p)));
            }
            else
            {
                // F^-1(p) = G^-1(1-p)
                return RationalApproximation(Math.Sqrt(-2.0 * Math.Log(1 - p)));
            }
        }



        static double cdf(double x)
        {
            double k = 1.0 / (1.0 + 0.2316419 * x);
            double k_sum = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + 1.330274429 * k))));

            if (x >= 0.0)
                return (1.0 - (1.0 / (Math.Pow(2 * Math.PI, 0.5))) * Math.Exp(-0.5 * x * x) * k_sum);
            else
                return 1.0 - cdf(-x);
        }

        // Inverse cumulative distribution function (aka the probit function)
        static double inv_cdf(double quantile)
        {
            // This is the Beasley-Springer-Moro algorithm which can 
            // be found in Glasserman [2004]. We won't go into the
            // details here, so have a look at the reference for more info
            double[] a = {2.50662823884,
                          -18.61500062529,
                          41.39119773534,
                         -25.44106049637};

            double[] b = {-8.47351093090,
                          23.08336743743,
                          -21.06224101826,
                          3.13082909833};

            double[] c = {0.3374754822726147,
                            0.9761690190917186,
                            0.1607979714918209,
                            0.0276438810333863,
                            0.0038405729373609,
                            0.0003951896511919,
                            0.0000321767881768,
                            0.0000002888167364,
                            0.0000003960315187};

            if (quantile >= 0.5 && quantile <= 0.92)
            {
                double num = 0.0;
                double denom = 1.0;

                for (int i = 0; i < 4; i++)
                {
                    num += a[i] * Math.Pow((quantile - 0.5), 2 * i + 1);
                    denom += b[i] * Math.Pow((quantile - 0.5), 2 * i);
                }
                return num / denom;

            }
            else if (quantile > 0.92 && quantile < 1)
            {
                double num = 0.0;

                for (int i = 0; i < 9; i++)
                {
                    num += c[i] * Math.Pow((Math.Log(-Math.Log(1 - quantile))), i);
                }
                return num;

            }
            else
                return -1.0 * inv_cdf(1 - quantile);
        }
    }

    public static class InverseCDF
    {
        public static double QNorm(double p, double mu, double sigma, bool lower_tail, bool log_p)
        {
            if (double.IsNaN(p) || double.IsNaN(mu) || double.IsNaN(sigma)) return (p + mu + sigma);
            double ans;
            bool isBoundaryCase = R_Q_P01_boundaries(p, double.NegativeInfinity, double.PositiveInfinity, lower_tail, log_p, out ans);
            if (isBoundaryCase) return (ans);
            if (sigma < 0) return (double.NaN);
            if (sigma == 0) return (mu);

            double p_ = R_DT_qIv(p, lower_tail, log_p);
            double q = p_ - 0.5;
            double r, val;

            if (Math.Abs(q) <= 0.425)  // 0.075 <= p <= 0.925
            {
                r = .180625 - q * q;
                val = q * (((((((r * 2509.0809287301226727 +
                           33430.575583588128105) * r + 67265.770927008700853) * r +
                         45921.953931549871457) * r + 13731.693765509461125) * r +
                       1971.5909503065514427) * r + 133.14166789178437745) * r +
                     3.387132872796366608)
                / (((((((r * 5226.495278852854561 +
                         28729.085735721942674) * r + 39307.89580009271061) * r +
                       21213.794301586595867) * r + 5394.1960214247511077) * r +
                     687.1870074920579083) * r + 42.313330701600911252) * r + 1.0);
            }
            else
            {
                r = q > 0 ? R_DT_CIv(p, lower_tail, log_p) : p_;
                r = Math.Sqrt(-((log_p && ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ? p : Math.Log(r)));

                if (r <= 5)              // <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11
                {
                    r -= 1.6;
                    val = (((((((r * 7.7454501427834140764e-4 +
                            .0227238449892691845833) * r + .24178072517745061177) *
                          r + 1.27045825245236838258) * r +
                         3.64784832476320460504) * r + 5.7694972214606914055) *
                       r + 4.6303378461565452959) * r +
                      1.42343711074968357734)
                     / (((((((r *
                              1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                             r + .0151986665636164571966) * r +
                            .14810397642748007459) * r + .68976733498510000455) *
                          r + 1.6763848301838038494) * r +
                         2.05319162663775882187) * r + 1.0);
                }
                else                     // very close to  0 or 1 
                {
                    r -= 5.0;
                    val = (((((((r * 2.01033439929228813265e-7 +
                            2.71155556874348757815e-5) * r +
                           .0012426609473880784386) * r + .026532189526576123093) *
                         r + .29656057182850489123) * r +
                        1.7848265399172913358) * r + 5.4637849111641143699) *
                      r + 6.6579046435011037772)
                     / (((((((r *
                              2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
                             r + 1.8463183175100546818e-5) * r +
                            7.868691311456132591e-4) * r + .0148753612908506148525)
                          * r + .13692988092273580531) * r +
                         .59983220655588793769) * r + 1.0);
                }
                if (q < 0.0) val = -val;
            }

            return (mu + sigma * val);
        }

        private static bool R_Q_P01_boundaries(double p, double _LEFT_, double _RIGHT_, bool lower_tail, bool log_p, out double ans)
        {
            if (log_p)
            {
                if (p > 0.0)
                {
                    ans = double.NaN;
                    return (true);
                }
                if (p == 0.0)
                {
                    ans = lower_tail ? _RIGHT_ : _LEFT_;
                    return (true);
                }
                if (p == double.NegativeInfinity)
                {
                    ans = lower_tail ? _LEFT_ : _RIGHT_;
                    return (true);
                }
            }
            else
            {
                if (p < 0.0 || p > 1.0)
                {
                    ans = double.NaN;
                    return (true);
                }
                if (p == 0.0)
                {
                    ans = lower_tail ? _LEFT_ : _RIGHT_;
                    return (true);
                }
                if (p == 1.0)
                {
                    ans = lower_tail ? _RIGHT_ : _LEFT_;
                    return (true);
                }
            }
            ans = double.NaN;
            return (false);
        }

        private static double R_DT_qIv(double p, bool lower_tail, bool log_p)
        {
            return (log_p ? (lower_tail ? Math.Exp(p) : -ExpM1(p)) : R_D_Lval(p, lower_tail));
        }

        private static double R_DT_CIv(double p, bool lower_tail, bool log_p)
        {
            return (log_p ? (lower_tail ? -ExpM1(p) : Math.Exp(p)) : R_D_Cval(p, lower_tail));
        }

        private static double R_D_Lval(double p, bool lower_tail)
        {
            return lower_tail ? p : 0.5 - p + 0.5;
        }

        private static double R_D_Cval(double p, bool lower_tail)
        {
            return lower_tail ? 0.5 - p + 0.5 : p;
        }
        private static double ExpM1(double x)
        {
            if (Math.Abs(x) < 1e-5)
                return x + 0.5 * x * x;
            else
                return Math.Exp(x) - 1.0;
        }
    }

    public class Test
    {
        public int sum = 0;

        public Test()
        {

        }

        public void ExecuteTest()
        {
            sum++;
        }
    }
        
}
