using System;
using OptimizationMethods.Optimization;
using OptimizationMethods.Optimization.SteepestDescentMethod;
using OptimizationMethods.Optimization.NonLinearConjugateGradient;
using OptimizationMethods.Optimization.LinearSystem;
using OptimizationMethods.Optimization.SequentialQuadraticProgramming;
using OptimizationMethods.Utils;

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
            Func<double[], double> testFunc1 = (x) =>
            {
                return Math.Pow(x[0], 4) - 8 * Math.Pow(x[0], 2) + 5;
            };

            Func<double[], double>[] dtestFunc1 = new Func<double[], double>[1];

            dtestFunc1[0] = (x) =>
            {
                return 4 * Math.Pow(x[0], 3) - 16 * x[0];
            };

            Func<double[], double> testConstr = (x) =>
            {
                return 5 - Math.Exp(x[0]) + 2.0 * Math.Pow(x[0] - 1, 2);
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

            Func<double[], double> bananaFunc = (x) =>
            {
                return Math.Pow(1 - x[0], 2) + 100 * Math.Pow(x[1] - x[0] * x[0], 2);
            };

            Func<double[], double> powell = (x) =>
            {
                return Math.Pow(x[0] + 10*x[1], 2) + 
                       5 * Math.Pow(x[2] - x[3], 2) +
                       Math.Pow(x[1] + 2 * x[2], 4) +
                       10 * Math.Pow(x[0] - x[3], 4);
            };

            OptVector[] ttt = new OptVector[3];
            ttt[0] = new OptVector(new double[3] { 1, 2 ,3 });
            ttt[1] = new OptVector(new double[3] { 4, 5, 6 });
            ttt[2] = new OptVector(new double[3] { 7, 8, 9 });

            var tr = OptVector.Transpose(ttt);

            //TestsSQP.Test0();
            TestsSQP.Test1();
            TestsSQP.Test2();
            TestsSQP.Test3();
            ////TestsSQP.Test4();
            ////TestsSQP.Test5();
            TestsSQP.Test6();
            TestsSQP.Test7();
            TestsSQP.Test8();
            TestsSQP.Test9();
            TestsSQP.Test10();
            TestsSQP.Test11();
            TestsSQP.Test12();
            TestsSQP.Test13();
            TestsSQP.Test14();
            TestsSQP.Test15();
            TestsSQP.Test16();
            TestsSQP.Test17();
            TestsSQP.Test18();
            TestsSQP.Test19();
            TestsSQP.Test20();
            TestsSQP.Test21();
            TestsSQP.Test22();
            TestsSQP.Test23();
            TestsSQP.Test24();
            TestsSQP.Test25();
            TestsSQP.Test26();
            TestsSQP.Test27();
            TestsSQP.Test28();
            //TestCGMethod();

            BFGS bfsg = new BFGS();

            var  result0 = bfsg.Solve(powell, new double[] { 3.0, -1.0, 0.0, 1.0 }, 10000);

            Console.WriteLine("Min " + powell(result0));

            NLCG gradient = new NLCG();

            var result = gradient.Solve(powell,   new double[] { 3.0, -1.0, 0.0, 1.0 }, 4000);

            Console.WriteLine("Min " + powell(result));

            SteepestDescent gradientSteep = new SteepestDescent();

            var result1 = gradientSteep.Solve(powell, new double[] { 3.0, -1.0, 0.0, 1.0 }, 10000);

            Console.WriteLine("Min " + powell(result1));

            //Variables result = SteepestDescent(testFunc, dTestFunc, 2, 20);

            for (int i = 0; i < result.Length; i++)
                Console.WriteLine("result " + result[i]);

            //Console.WriteLine("ver " + testFunc(result));

            Console.ReadLine();
        }


        
        static void TestFunc()
        {
            OptVector a = new OptVector(new double[] { 1, 2 });
            OptVector b = new OptVector(new double[] { 3, 4 });

            OptVector c = new OptVector(new double[] { 1, 3 });
            OptVector d = new OptVector(new double[] { 4, 7 });

            var res = OptVector.Mult(a, b);

            var res1 = OptVector.SubtractFromIdentity(res);

            var res2 = OptVector.Div(res, 2);

            OptVector[] aa = new OptVector[] { a, b };
            OptVector[] bb = new OptVector[] { c, d };

            var res3 = OptVector.Mult(res, bb);
            var res4 = OptVector.Sum(res, bb);

            var res5 = OptVector.Mult(res, c);

            
        }

        

        static void TestCGMethod()
        {
            MINRES minres = new MINRES();
            CGMethod cg = new CGMethod();

            OptVector[] A2 = new OptVector[3];
            A2[0] = new OptVector(new double[] { 2, 1, 3 });
            A2[1] = new OptVector(new double[] { 2, 6, 8 });
            A2[2] = new OptVector(new double[] { 6, 8, 18 });

            OptVector b2 = new OptVector(new double[] { 1, 3, 5 });

            var minresSol = minres.Solve(A2, b2, new OptVector(new double[3]), 50);
            var sol2 = cg.Solve(A2, b2, new OptVector(new double[3]), 50);

            OptVector[] A = new OptVector[3];
            A[0] = new OptVector(new double[] { 2, 0, 1 });
            A[1] = new OptVector(new double[] { 1, 6, 0 });
            A[2] = new OptVector(new double[] { 3, 2, 3 });

            OptVector x = new OptVector(new double[] { 2, 5, 7 });
            
            var sol = cg.Solve(A, x, new OptVector(new double[3]), 50);
            var solm = minres.Solve(A, x, new OptVector(new double[3]), 50);

            OptVector[] A1 = new OptVector[3];
            A1[0] = new OptVector(new double[] { 3, 1, -6 });
            A1[1] = new OptVector(new double[] { 2, 1, -5 });
            A1[2] = new OptVector(new double[] { 6, -3, 3 });

            OptVector b1 = new OptVector(new double[] { -10, -8, 0 });

            var sol1 = cg.Solve(A1, b1, new OptVector(new double[3]), 200);
            var solm1 = minres.Solve(A1, b1, new OptVector(new double[3]), 200);

            OptVector diff = b1 - OptVector.Mult(A1, solm1);
            OptVector diff1 = b1 - OptVector.Mult(A1, sol1);
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
