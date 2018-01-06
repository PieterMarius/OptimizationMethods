using ConsoleApplication1.Optimization;
using ConsoleApplication1.Optimization.SequentialQuadraticProgramming;
using System;
using System.Collections.Generic;

namespace ConsoleApplication1
{
    public static class TestsSQP
    {
        static void Test0()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {
                return Math.Pow(x[0], 2) +
                       Math.Pow(x[1], 2);
            };

            List<Func<double[], double>> eqConstraint = new List<Func<double[], double>>();

            Func<double[], double> eqConstraint1 = (x) =>
            {
                return -x[0] -
                       x[1] + 2;
            };

            eqConstraint.Add(eqConstraint1);

            var res = quadraticProgramming.Minimize(f, eqConstraint, null, new double[] { 0, 0 }, 50);
        }

        public static void Test1()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {
                return Math.Pow(x[0] - 1, 2) +
                       Math.Pow(x[1] - 3, 2);
            };

            List<Func<double[], double>> eqConstraint = new List<Func<double[], double>>();

            Func<double[], double> eqConstraint1 = (x) =>
            {
                return -x[0] +
                       Math.Pow(x[1], 2) - 1;
            };

            eqConstraint.Add(eqConstraint1);

            var res = quadraticProgramming.Minimize(f, eqConstraint, null, new double[] { 0, 0 }, 100);
        }

        /// <summary>
        /// x1 = 5/3, x2 = 1/3
        /// </summary>
        public static void Test2()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {
                return -(5 - Math.Pow(x[0] - 2, 2) -
                       2 * Math.Pow(x[1] - 1, 2));
            };

            List<Func<double[], double>> eqConstraint = new List<Func<double[], double>>();

            Func<double[], double> eqConstraint1 = (x) =>
            {
                return -x[0] -
                       4 * x[1] + 3;
            };

            eqConstraint.Add(eqConstraint1);

            var res = quadraticProgramming.Minimize(f, eqConstraint, null, new double[] { 0, 0 }, 100);

            double viol = eqConstraint1(res);
            double viol1 = eqConstraint1(new double[] { 5.0 / 3.0, 1.0 / 3.0 });
        }

        /// <summary>
        /// x1= 1.0, x2 = 0.0, x3 = 0.0
        /// </summary>
        public static void Test3()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {
                return x[0] + x[1] + Math.Pow(x[2], 2);
            };

            List<Func<double[], double>> eqConstraint = new List<Func<double[], double>>();

            Func<double[], double> eqConstraint1 = (x) =>
            {
                return x[0] - 1;
            };

            Func<double[], double> eqConstraint2 = (x) =>
            {
                return Math.Pow(x[0], 2) + Math.Pow(x[1], 2) - 1;
            };

            eqConstraint.Add(eqConstraint1);
            eqConstraint.Add(eqConstraint2);

            var res = quadraticProgramming.Minimize(f, eqConstraint, null, new double[] { 0.1, 0.1, 0.1 }, 1000);
        }

        public static void Test4()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {

                return 400 * Math.Pow(x[0], 2) +
                       800 * Math.Pow(x[1], 2) +
                       200 * x[0] * x[1] +
                       1600 * Math.Pow(x[2], 2) +
                       400 * x[1] * x[2];
            };

            List<Func<double[], double>> eqConstraint = new List<Func<double[], double>>();

            Func<double[], double> eqConstraint1 = (x) =>
            {
                return 1.2 - x[0] - x[1] - 1.5 * x[2];
            };

            Func<double[], double> eqConstraint2 = (x) =>
            {
                return 1.0 - x[0] - x[1] - x[2];
            };

            eqConstraint.Add(eqConstraint2);
            eqConstraint.Add(eqConstraint1);

            var res = quadraticProgramming.Minimize(f, eqConstraint, null, new double[] { 0, 0, 0 }, 50);
        }

        public static void Test5()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {

                return Math.Pow(x[0] - 2, 2) +
                       2 * Math.Pow(x[1] - 1, 2);
            };

            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();

            Func<double[], double> inqConstraint1 = (x) =>
            {
                return 3 - x[0] - 4 * x[1];
            };

            Func<double[], double> inqConstraint2 = (x) =>
            {
                return x[0] - x[1];
            };

            inqConstraint.Add(inqConstraint1);
            inqConstraint.Add(inqConstraint2);

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { 1, 1 }, 50);
        }

        /// <summary>
        /// x1 = 8.5, x2 = 8.75, x3 = 17.25
        /// </summary>
        public static void Test6()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {

                return -(x[0] * (30 - x[0]) + x[1] * (50 - 2 * x[1]) - 3 * x[0] - 5 * x[1] - 10 * x[2]);
            };

            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();

            Func<double[], double> inqConstraint1 = (x) =>
            {
                return x[0] + x[1] - x[2];
            };

            Func<double[], double> inqConstraint2 = (x) =>
            {
                return x[2] - 17.25;
            };

            inqConstraint.Add(inqConstraint1);
            inqConstraint.Add(inqConstraint2);

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { 0, 0, 0 }, 200);
        }

        /// <summary>
        /// x1 = 0.688, x2 = 0.883
        /// </summary>
        public static void Test7()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {

                return -(Math.Sin(x[0]) * Math.Cos(x[1]) + Math.Cos(x[0]) * Math.Sin(x[1]));
            };

            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();

            Func<double[], double> inqConstraint1 = (x) =>
            {
                return -x[0] - x[1];
            };

            Func<double[], double> inqConstraint2 = (x) =>
            {
                return x[0] + x[1] - Math.PI;
            };

            inqConstraint.Add(inqConstraint1);
            inqConstraint.Add(inqConstraint2);

            List<Func<double[], double>> eqConstraint = new List<Func<double[], double>>();

            Func<double[], double> eqConstraint1 = (x) =>
            {
                return x[0] - Math.Pow(x[1], 3);
            };

            eqConstraint.Add(eqConstraint1);

            var res = quadraticProgramming.Minimize(f, eqConstraint, inqConstraint, new double[] { 0.0, 0.0 }, 1000);
        }

        /// <summary>
        /// x1 = -1.0, x2 = -1.0
        /// </summary>
        public static void Test8()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {

                return x[0] + x[1];
            };

            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();

            Func<double[], double> inqConstraint1 = (x) =>
            {
                return -2 + Math.Pow(x[0], 2) + Math.Pow(x[1], 2);
            };

            inqConstraint.Add(inqConstraint1);

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { 0, 0 }, 1000);
        }


        /// <summary>
        /// Rosenbrock's function (x1 = 0.4149, x2 = 0.1701)
        /// </summary>
        public static void Test9()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {

                return Math.Pow(1 - x[0], 2) + 100 * Math.Pow(x[1] - Math.Pow(x[0], 2), 2);
            };

            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();

            Func<double[], double> inqConstraint1 = (x) =>
            {
                return x[0] + 2 * x[1] - 1;
            };

            inqConstraint.Add(inqConstraint1);

            List<Func<double[], double>> eqConstraint = new List<Func<double[], double>>();

            Func<double[], double> eqConstraint1 = (x) =>
            {
                return 2 * x[0] + x[1] - 1;
            };

            eqConstraint.Add(eqConstraint1);

            var res = quadraticProgramming.Minimize(f, eqConstraint, inqConstraint, new double[] { 0.5, 0 }, 100);
        }

        /// <summary>
        /// Bound constraint (x1 = 1, x2 = 2)
        /// </summary>
        public static void Test10()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {

                return 1 + x[0] / (1 + x[1]) - 3 * x[0] * x[1] + x[1] * (1 + x[0]);
            };

            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();

            Func<double[], double> inqConstraint1 = (x) =>
            {
                return x[0] - 1;
            };

            Func<double[], double> inqConstraint2 = (x) =>
            {
                return x[1] - 2;
            };

            Func<double[], double> inqConstraint3 = (x) =>
            {
                return -x[0];
            };


            Func<double[], double> inqConstraint4 = (x) =>
            {
                return -x[1];
            };

            inqConstraint.Add(inqConstraint1);
            inqConstraint.Add(inqConstraint2);
            inqConstraint.Add(inqConstraint3);
            inqConstraint.Add(inqConstraint4);


            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { 0.5, 1.0 }, 100);
        }


        /// <summary>
        /// Rosenbrock's function non linear constraints (x1 = 0.5, x2 = 0.25)
        /// </summary>
        public static void Test11()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {

                return Math.Pow(1.0 - x[0], 2) + 100.0 * Math.Pow(x[1] - Math.Pow(x[0], 2), 2);
            };

            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();

            Func<double[], double> eqConstraint1 = (x) =>
            {
                return Math.Pow(x[0] - (1.0 / 3.0), 2) + Math.Pow(x[1] - (1.0 / 3.0), 2) - Math.Pow(1.0 / 3.0, 2);
            };

            //Func<Vector, double> inqConstraint1 = (x) =>
            //{
            //    return x[0] - 0.5;
            //};

            //Func<Vector, double> inqConstraint2 = (x) =>
            //{
            //    return x[1] - 0.8;
            //};

            //Func<Vector, double> inqConstraint3 = (x) =>
            //{
            //    return -x[0];
            //};

            //Func<Vector, double> inqConstraint4 = (x) =>
            //{
            //    return -x[1] + 0.2;
            //};

            inqConstraint.Add(eqConstraint1);
            //inqConstraint.Add(inqConstraint1);
            //inqConstraint.Add(inqConstraint2);
            //inqConstraint.Add(inqConstraint3);
            //inqConstraint.Add(inqConstraint4);


            double[] upperBound = new double[] { 0.5, 0.8 };
            double[] lowerBound = new double[] { 0.0, 0.2 };

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, lowerBound, upperBound, new double[] { 0.25, 0.25 }, 1000);
        }


        /// <summary>
        /// x1 = 1.4, x2 = 1.7
        /// </summary>
        public static void Test12()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {

                return Math.Pow(x[0] - 1.0, 2) + Math.Pow(x[1] - 2.5, 2);
            };

            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();

            Func<double[], double> inqConstraint1 = (x) =>
            {
                return -x[0] + 2 * x[1] - 2.0;
            };

            Func<double[], double> inqConstraint2 = (x) =>
            {
                return x[0] + x[1] - 6.0;
            };

            Func<double[], double> inqConstraint3 = (x) =>
            {
                return x[0] - 2.0 * x[1] - 2.0;
            };

            Func<double[], double> inqConstraint4 = (x) =>
            {
                return -x[0];
            };

            Func<double[], double> inqConstraint5 = (x) =>
            {
                return -x[1];
            };

            inqConstraint.Add(inqConstraint1);
            inqConstraint.Add(inqConstraint2);
            inqConstraint.Add(inqConstraint3);
            inqConstraint.Add(inqConstraint4);
            inqConstraint.Add(inqConstraint5);


            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { 2.0, 0.0 }, 100);
        }

        /// <summary>
        /// Rosenbrock's function non linear constraints (x1 = 0.5, x2 = 0.25)
        /// </summary>
        public static void Test13()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {

                return Math.Pow(1.0 - x[0], 2) + 100.0 * Math.Pow(x[1] - Math.Pow(x[0], 2), 2);
            };

            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();

            Func<double[], double> inqConstraint1 = (x) =>
            {
                return x[0] - 0.5;
            };

            Func<double[], double> inqConstraint2 = (x) =>
            {
                return x[1] - 0.8;
            };

            Func<double[], double> inqConstraint3 = (x) =>
            {
                return -x[0];
            };

            Func<double[], double> inqConstraint4 = (x) =>
            {
                return -x[1] + 0.2;
            };

            Func<double[], double> eqConstraint1 = (x) =>
            {
                return Math.Pow(x[0] - (1.0 / 3.0), 2) + Math.Pow(x[1] - (1.0 / 3.0), 2) - Math.Pow(1.0 / 3.0, 2);
            };

            inqConstraint.Add(inqConstraint1);
            inqConstraint.Add(inqConstraint2);
            inqConstraint.Add(eqConstraint1);
            inqConstraint.Add(inqConstraint3);
            inqConstraint.Add(inqConstraint4);

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { 0.25, 0.25 }, 1000);
        }

        /// <summary>
        /// Rosenbrock's function non linear constraints (x1 = 0.7864, x2 = 0.6177)
        /// </summary>
        public static void Test14()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {

                return Math.Pow(1.0 - x[0], 2) + 100.0 * Math.Pow(x[1] - Math.Pow(x[0], 2), 2);
            };

            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();

            Func<double[], double> eqConstraint1 = (x) =>
            {
                return Math.Pow(x[0], 2) + Math.Pow(x[1], 2) - 1;
            };

            inqConstraint.Add(eqConstraint1);

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { 0.0, 0.0 }, 1000);

            double min = f(res);
            double eqc = eqConstraint1(res);

            double min1 = f(new double[] { 0.7864, 0.6177 });
            double eqc2 = eqConstraint1(new double[] { 0.7864, 0.6177 });
        }

        public static void Test15()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {

                return -(2 * x[0] * x[1] + 2 * x[0] - Math.Pow(x[0], 2) - 2 * Math.Pow(x[1], 2));
            };

            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();

            Func<double[], double> inqConstraint1 = (x) =>
            {
                return -x[1] + 1.0;
            };

            inqConstraint.Add(inqConstraint1);

            List<Func<double[], double>> eqConstraint = new List<Func<double[], double>>();

            Func<double[], double> eqConstraint1 = (x) =>
            {
                return Math.Pow(x[0], 3) - x[1];
            };

            eqConstraint.Add(eqConstraint1);

            var res = quadraticProgramming.Minimize(f, eqConstraint, inqConstraint, new double[] { 0.0, 0.0 }, 1000);
        }


        /// <summary>
        /// x1 = 2.0, x2 = 0.0
        /// </summary>
        public static void Test16()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {
                return 0.01 * Math.Pow(x[0], 2) + Math.Pow(x[1], 2) - 100.0;
            };

            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();

            Func<double[], double> inqConstraint1 = (x) =>
            {
                return x[0] - 50;
            };

            Func<double[], double> inqConstraint2 = (x) =>
            {
                return x[1] - 50;
            };

            Func<double[], double> inqConstraint3 = (x) =>
            {
                return -x[0] + 2.0;
            };

            Func<double[], double> inqConstraint4 = (x) =>
            {
                return -x[1] - 50;
            };

            Func<double[], double> inqConstraint5 = (x) =>
            {
                return -10 * x[0] + x[1] + 10.0;
            };

            inqConstraint.Add(inqConstraint1);
            inqConstraint.Add(inqConstraint2);
            inqConstraint.Add(inqConstraint3);
            inqConstraint.Add(inqConstraint4);
            inqConstraint.Add(inqConstraint5);

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { -1.0, -1.0 }, 100);
        }

        /// <summary>
        /// Rosenbrock's function non linear constraints (x1 = 0.9072, x2 = 0.8228)
        /// </summary>
        public static void Test17()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {

                return Math.Pow(1.0 - x[0], 2) + 100.0 * Math.Pow(x[1] - Math.Pow(x[0], 2), 2);
            };

            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();

            Func<double[], double> eqConstraint1 = (x) =>
            {
                return Math.Pow(x[0], 2) + Math.Pow(x[1], 2) - 1.5;
            };

            inqConstraint.Add(eqConstraint1);

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { 1.0, -0.5 }, 200);

            double test = eqConstraint1(res);
            double min = f(res);

            double test1 = eqConstraint1(new double[] { 0.9072, 0.8228 });
            double min1 = f(new double[] { 0.9072, 0.8228 });
        }

        public static void Test18()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {

                return Math.Pow(x[0], 2) + Math.Pow(x[1], 2) + Math.Pow(x[2], 2);
            };
            //−x1 − x2 exp (x3y) − exp(2y) + 2 exp(4y) ≥ 0 for all y ∈ [0, 1]

            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();


            int nConstraints = 100;
            double step = 1.0 / nConstraints;
            double aStep = 0.0;

            for (int i = 0; i < nConstraints; i++)
            {
                aStep += step;

                double inner = -aStep;

                Func<double[], double> eqConstraint1 = (x) =>
                {
                    return x[0] + x[1] * Math.Exp(x[2] * inner) + Math.Exp(2 * inner) - 2 * Math.Exp(4 * inner);
                };

                inqConstraint.Add(eqConstraint1);
                //bool tt = false;
                //if (eqConstraint1(new double[] { 1.0, -1.0, 2.0 }) > 0.0)
                //    tt = true;
            }

            double min = f(new double[] { 1.0, -1.0, 2.0 });

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { 1.0, -1.0, 2.0 }, 30000);

            double min1 = f(res);

            foreach (var func in inqConstraint)
            {
                double test = func(res);
                var t = false;
                if (test > 0)
                    t = true;
            }
        }

        /// <summary>
        /// x1 = 1.6667, x2 = 2.11111
        /// </summary>
        public static void Test19()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {

                return Math.Pow(x[0], 4) + Math.Pow(x[1], 4);
            };


            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();

            Func<double[], double> eqConstraint1 = (x) =>
            {
                return Math.Pow(x[0], 2) - x[0] - x[1] + 1;
            };

            Func<double[], double> eqConstraint2 = (x) =>
            {
                return Math.Pow(x[0], 2) - 4 * x[0] - x[1] + 6;
            };

            Func<double[], double> eqConstraint3 = (x) =>
            {
                return Math.Pow(x[0], 2) - 3 * x[0] + x[1] - 2;
            };

            inqConstraint.Add(eqConstraint1);
            inqConstraint.Add(eqConstraint2);
            inqConstraint.Add(eqConstraint3);

            foreach (var func in inqConstraint)
            {
                double test = func(new double[] { 1.66667, 2.11111 });
            }

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { -4.0, -4.0 }, 10000);
        }

        public static void Test20()
        {
            SQP quadraticProgramming = new SQP();

            int nVar = 16;

            Func<double[], double> f = (x) =>
            {
                double sum = 0.0;
                for (int i = 0; i < nVar - 1; i++)
                {
                    sum += 100.0 * Math.Pow(x[i + 1] - x[i] * x[i], 2) + Math.Pow(x[i] - 1.0, 2);
                }
                return 0.5 * sum;
            };

            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();

            Func<double[], double> inqConstr = (x) =>
            {
                double sum = 0.0;
                for (int i = 0; i < nVar - 1; i++)
                {
                    sum += 1.1 - Math.Pow(x[i] - 2.0, 3) - x[i + 1];
                }
                return -sum;
            };

            inqConstraint.Add(inqConstr);

            double[] startValue = new double[nVar];
            startValue[0] = 4.0;
            for (int i = 1; i < nVar; i++)
            {
                startValue[i] = 4.0;
            }

            double testStart = inqConstr(startValue);

            //var solver = new BFGS(1E-50, true);
            //var solbf = solver.Solve(f, startValue, 5000);
            //var testbf = f(solbf);

            var stmin = f(startValue);
            var res = quadraticProgramming.Minimize(f, null, inqConstraint,startValue, 1000);
            //var res = quadraticProgramming.Minimize(f, null, null, startValue, 10000);

            double test = inqConstr(res);
            double min = f(res);
        }

        public static void Test21()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {

                return Math.Pow(x[0] - 1.0, 2) + Math.Pow(x[0] - x[1], 2) +
                       Math.Pow(x[2] - 1.0, 2) + Math.Pow(x[3] - 1.0, 4) +
                       Math.Pow(x[4] - 1.0, 6);
            };


            List<Func<double[], double>> eqConstraint = new List<Func<double[], double>>();

            Func<double[], double> eqConstraint1 = (x) =>
            {
                return Math.Pow(x[0], 2) * x[3] + Math.Sin(x[3] - x[4]) - 2 * Math.Sqrt(2);
            };

            Func<double[], double> eqConstraint2 = (x) =>
            {
                return x[1] + Math.Pow(x[2], 4) * Math.Pow(x[3], 2) - 8 - Math.Sqrt(2);
            };

            eqConstraint.Add(eqConstraint1);
            eqConstraint.Add(eqConstraint2);

            double[] startValue = new double[5];
            for (int i = 0; i < 5; i++)
            {
                startValue[i] = 2.0;
            }
            

            var res = quadraticProgramming.Minimize(f, eqConstraint, null, startValue, 10000);
        }

        public static void Test22()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {

                return (x[0] * x[0]) + (x[1] * x[1]) + x[0] * x[1] - 14.0 * x[0] - 16 * x[1] +
                        Math.Pow(x[2] - 10.0, 2) + 4 * Math.Pow(x[3] - 5, 2) + Math.Pow(x[4] - 3.0, 2) +
                       2 * Math.Pow(x[5] - 1.0, 2) + 5 * x[6] * x[6] + 7 * Math.Pow(x[7] - 11.0, 2) +
                       2 * Math.Pow(x[8] - 10.0, 2) + Math.Pow(x[9] - 7.0, 2) + 45.0;
            };


            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();

            Func<double[], double> inqConstraint1 = (x) =>
            {
                return -(105.0 - 4 * x[0] - 5 * x[1] + 3.0 * x[7] - 9.0 * x[7]);
            };

            Func<double[], double> inqConstraint2 = (x) =>
            {
                return -(- 10 * x[0] + 8 * x[1] + 17.0 * x[6] - 2.0 * x[7]);
            };

            Func<double[], double> inqConstraint3 = (x) =>
            {
                return -(8 * x[0] - 2 * x[1] - 5.0 * x[8] + 2.0 * x[9] +12);
            };

            Func<double[], double> inqConstraint4 = (x) =>
            {
                return -(-3.0 * Math.Pow(x[0] - 2.0, 2) - 4.0 * Math.Pow(x[1] - 3.0, 2) - 2.0 * x[2] * x[2] + 7 * x[3] + 120);
            };

            Func<double[], double> inqConstraint5 = (x) =>
            {
                return -(-5.0 * x[0]*x[0] - 8.0 * x[1] -Math.Pow(x[2]-6.0,2) + 2.0 * x[3]+ 40);
            };

            Func<double[], double> inqConstraint6 = (x) =>
            {
                return -(0.5 * Math.Pow(x[0] - 8.0, 2) - 2.0 * Math.Pow(x[1] - 4.0, 2) - 3.0 * x[4] * x[4] + x[5] + 30.0);
            };

            Func<double[], double> inqConstraint7 = (x) =>
            {
                return -(-x[0] * x[0] - 2.0 * Math.Pow(x[1] - 2.0, 2) + 2.0 * x[0] * x[1] - 14.0 * x[4] + 6.0 * x[5]);
            };

            Func<double[], double> inqConstraint8 = (x) =>
            {
                return -(3 * x[0] - 6.0 * x[1] - 12.0 * Math.Pow(x[8] - 8.0, 2) + 7 * x[9]);
            };

            inqConstraint.Add(inqConstraint1);
            inqConstraint.Add(inqConstraint2);
            inqConstraint.Add(inqConstraint3);
            inqConstraint.Add(inqConstraint4);
            inqConstraint.Add(inqConstraint5);
            inqConstraint.Add(inqConstraint6);
            inqConstraint.Add(inqConstraint7);
            inqConstraint.Add(inqConstraint8);

            double[] startValue = new double[10] { 2.0, 3.0, 5.0, 5.0, 1.0, 2.0, 7.0, 3.0, 6.0, 10.0 };
            double fval = f(startValue);
            
            var res = quadraticProgramming.Minimize(f, null, inqConstraint, startValue, 10000);

            double fval1 = f(res);

            foreach(var item in inqConstraint)
            {
                bool test = false;
                double vv = item(res);
                if (vv > 0)
                    test = true;
            }
        }
    }
}
