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

            var res = quadraticProgramming.Minimize(f, eqConstraint, null, new double[] { 0.1, 0.1, 0.1 }, 3000);
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

            var res = quadraticProgramming.Minimize(f, eqConstraint, inqConstraint, new double[] { 0.1, 0.1 }, 1000);
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
            SQP quadraticProgramming = new SQP(1E-30, true);

            Func<double[], double> f = (x) =>
            {

                return 1.0 + x[0] / (1.0 + x[1]) - 3.0 * x[0] * x[1] + x[1] * (1.0 + x[0]);
            };

            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();

            Func<double[], double> inqConstraint1 = (x) =>
            {
                return x[0] - 1.0;
            };

            Func<double[], double> inqConstraint2 = (x) =>
            {
                return x[1] - 2.0;
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

            var min = f(new double[] { 1, 2 });

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { 0.5, 1.0 }, 1000);

            var min1 = f(res);
        }


        /// <summary>
        /// Rosenbrock's function non linear constraints (x1 = 0.5, x2 = 0.25)
        /// </summary>
        public static void Test11()
        {
            SQP quadraticProgramming = new SQP(1E-25,true);

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
                return -x[0] + 2.0 * x[1] - 2.0;
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
                return Math.Pow(x[0], 2) + Math.Pow(x[1], 2) - 1.0;
            };

            inqConstraint.Add(eqConstraint1);

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { 0.1, 0.1 }, 5000);

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

            var res = quadraticProgramming.Minimize(f, eqConstraint, inqConstraint, new double[] { 0.0, 0.0 }, 3000);
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

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { 3.0, -1.0 }, 1000);
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

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { 1.0, -0.5 }, 1000);

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

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { 1.0, -1.0, 2.0 }, 3000);

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

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { -4.0, -4.0 }, 2000);
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
        //n 77
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
            

            var res = quadraticProgramming.Minimize(f, eqConstraint, null, startValue, 1000);

            double min = f(res);

            foreach (var func in eqConstraint)
            {
                double test = func(res);
                var t = false;
                if (test > 0)
                    t = true;
            }


        }

        //Test 113
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
                return -(105.0 - 4.0 * x[0] - 5.0 * x[1] + 3.0 * x[6] - 9.0 * x[7]);
            };

            Func<double[], double> inqConstraint2 = (x) =>
            {
                return -(- 10.0 * x[0] + 8.0 * x[1] + 17.0 * x[6] - 2.0 * x[7]);
            };

            Func<double[], double> inqConstraint3 = (x) =>
            {
                return -(8.0 * x[0] - 2.0 * x[1] - 5.0 * x[8] + 2.0 * x[9] +12.0);
            };

            Func<double[], double> inqConstraint4 = (x) =>
            {
                return -(-3.0 * Math.Pow(x[0] - 2.0, 2) - 4.0 * Math.Pow(x[1] - 3.0, 2) - 2.0 * x[2] * x[2] + 7.0 * x[3] + 120.0);
            };

            Func<double[], double> inqConstraint5 = (x) =>
            {
                return -(-5.0 * x[0] * x[0] - 8.0 * x[1] - Math.Pow(x[2] - 6.0, 2) + 2.0 * x[3] + 40.0);
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
                return -(3.0 * x[0] - 6.0 * x[1] - 12.0 * Math.Pow(x[8] - 8.0, 2) + 7.0 * x[9]);
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
            //double[] startValue = new double[10] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            double fval = f(startValue);
            
            var res = quadraticProgramming.Minimize(f, null, inqConstraint, startValue, 15000);

            double fval1 = f(res);

            foreach(var item in inqConstraint)
            {
                bool test = false;
                double vv = item(res);
                if (vv > 0)
                    test = true;
            }
        }

        //n 106
        public static void Test23()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {

                return (x[0] + x[1] + x[2]);
            };


            List<Func<double[], double>> inqConstraint = new List<Func<double[], double>>();

            Func<double[], double> inqConstraint1 = (x) =>
            {
                return -(1.0 - 0.0025 * (x[3] + x[5]));
            };

            Func<double[], double> inqConstraint2 = (x) =>
            {
                return -(1.0 - 0.0025 * (x[4] + x[6] - x[3]));
            };

            Func<double[], double> inqConstraint3 = (x) =>
            {
                return -(1.0 - 0.01 * (x[7] - x[4]));
            };

            Func<double[], double> inqConstraint4 = (x) =>
            {
                return -(x[0] * x[5] - 833.33252 * x[3] - 100.0 * x[0] + 83333.333);
            };

            Func<double[], double> inqConstraint5 = (x) =>
            {
                return -(x[1] * x[6] - 1250.0 * x[4] - x[1] * x[3] + 1250.0 * x[3]);
            };

            Func<double[], double> inqConstraint6 = (x) =>
            {
                return -(x[2] * x[7] - 1250000 - x[2] * x[4] + 2500.0 * x[4]);
            };

            inqConstraint.Add(inqConstraint1);
            inqConstraint.Add(inqConstraint2);
            inqConstraint.Add(inqConstraint3);
            inqConstraint.Add(inqConstraint4);
            inqConstraint.Add(inqConstraint5);
            inqConstraint.Add(inqConstraint6);

            double[] upperBound = new double[] { 10000.0, 10000.0, 10000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0 };
            double[] lowerBound = new double[] { 100.0, 1000.0, 1000.0, 10.0, 10.0, 10.0, 10.0, 10.0 };

            var boundsConstraints = CreateBoundsConstraints(lowerBound, upperBound);
            inqConstraint.AddRange(boundsConstraints);
            //double[] initalGuess = GetInitialGuess(inqConstraint, null, lowerBound, upperBound, 8);

            //foreach (var item in inqConstraint)
            //{
            //    bool test = false;
            //    double vv = item(initalGuess);
            //    if (vv > 0)
            //        test = true;
            //}
            
            //double[] startValue = new double[] { 5000.0, 5000.0, 5000.0, 200.0, 350.0, 150.0, 225.0, 425.0 };
            //StartSolution(inqConstraint,null, startValue);
            double[] startValue = new double[] { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
            double fval = f(startValue);

            

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, lowerBound, upperBound, startValue, 6000);

            double fval1 = f(res);

            foreach (var item in inqConstraint)
            {
                bool test = false;
                double vv = item(res);
                if (vv > 0)
                    test = true;
            }
        }

        //n 110
        public static void Test24()
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {
                double sum = 0.0;
                double prod = 1.0;
                for (int i = 0; i < 10; i++)
                {
                    sum += Math.Pow(Math.Log(x[i] - 2.0), 2) + Math.Pow(Math.Log(10.0 - x[i]), 2);
                    prod = prod * x[i];
                }

                sum -= Math.Pow(prod, 0.2);

                
                return sum;
            };


            double[] upperBound = new double[10]; 
            double[] lowerBound = new double[10];
            double[] startValue = new double[10];

            for (int i = 0; i < lowerBound.Length; i++)
            {
                lowerBound[i] = 2.001;
                upperBound[i] = 9.9;
                startValue[i] = 9.0;
            }

                       
            double fval = f(startValue);

            
            var res = quadraticProgramming.Minimize(f, null, null, lowerBound, upperBound, startValue, 2000);

            double fval1 = f(res);

            
        }

        //Dixon and price(48)
        public static void Test25()
        {
            SQP quadraticProgramming = new SQP();

            int nvalue = 4;

            Func<double[], double> f = (x) =>
            {
                double sum = Math.Pow(x[0] - 1.0, 2);

                for (int i = 1; i < nvalue - 1; i++)
                {
                    sum += (i + 1) * Math.Pow(2.0 * x[i] * x[i] - x[i - 1], 2);
                }
                               
                return sum;
            };


            double[] upperBound = new double[nvalue];
            double[] lowerBound = new double[nvalue];
            double[] startValue = new double[nvalue];

            for (int i = 0; i < lowerBound.Length; i++)
            {
                lowerBound[i] = -10.0;
                upperBound[i] = 10.0;
                startValue[i] = 0.1;
            }

            double fval = f(startValue);
            
            var res = quadraticProgramming.Minimize(f, null, null, lowerBound, upperBound, startValue, 2000);

            double fval1 = f(res);


        }

        //Jennrich function(67)
        public static void Test26()
        {
            SQP quadraticProgramming = new SQP(1E-50, true);

            int nvalue = 10;

            Func<double[], double> f = (x) =>
            {
                double sum = 0.0;

                for (int i = 1; i < nvalue+1; i++)
                {
                    sum += Math.Pow(2.0 + 2.0 * i - (Math.Exp(i * x[0]) + Math.Exp(i * x[1])), 2);
                }

                return sum;
            };


            double[] upperBound = new double[2];
            double[] lowerBound = new double[2];
            double[] startValue = new double[2];

            for (int i = 0; i < lowerBound.Length; i++)
            {
                lowerBound[i] = -1.0;
                upperBound[i] = 1.0;
                startValue[i] = -1.0;
            }

            double fval = f(startValue);

            var res = quadraticProgramming.Minimize(f, null, null, lowerBound, upperBound, startValue, 2000);

            double fval1 = f(res);


        }
        //TangFunction 143
        public static void Test27()
        {
            SQP quadraticProgramming = new SQP(1E-50,true);

            int nvalue = 2;

            Func<double[], double> f = (x) =>
            {
                double sum = 0.0;

                for (int i = 0; i < nvalue; i++)
                {
                    sum += Math.Pow(x[i], 4) - 16.0 * x[i] * x[i] + 5.0 * x[i];
                }

                return 0.5 * sum;
            };


            double[] upperBound = new double[nvalue];
            double[] lowerBound = new double[nvalue];
            double[] startValue = new double[nvalue];

            for (int i = 0; i < lowerBound.Length; i++)
            {
                lowerBound[i] = -5.0;
                upperBound[i] = 5.0;
                startValue[i] = 5.0;
            }

            double fval = f(startValue);

            var res = quadraticProgramming.Minimize(f, null, null, lowerBound, upperBound, startValue, 7000);

            double fval1 = f(res);


        }

        //119
        public static void Test28()
        {
            SQP quadraticProgramming = new SQP();

            int nvalue = 16;

            double[,] m = new double[nvalue, nvalue];
            for (int i = 0; i < nvalue; i++)
                m[i, i] = 1.0;

            m[0, 3] = 1.0; m[0, 6] = 1.0; m[0, 7] = 1.0; m[0, 15] = 1.0;

            m[1, 2] = 1.0; m[1, 6] = 1.0; m[1, 9] = 1.0;

            m[2, 6] = 1.0; m[2, 8] = 1.0; m[2, 9] = 1.0; m[2, 13] = 1.0;

            m[3, 6] = 1.0; m[3, 10] = 1.0; m[3, 14] = 1.0;

            m[4, 5] = 1.0; m[4, 9] = 1.0; m[4, 11] = 1.0; m[4, 15] = 1.0;

            m[5, 7] = 1.0; m[5, 14] = 1.0;

            m[6, 10] = 1.0;  m[6, 12] = 1.0;

            m[7, 9] = 1.0; m[7, 14] = 1.0;

            m[8, 11] = 1.0; m[8, 15] = 1.0;

            m[9, 13] = 1.0;
            
            m[10, 12] = 1.0;

            m[11, 13] = 1.0;

            m[12, 13] = 1.0;

            Func<double[], double> f = (x) =>
            {
                double sum = 0.0;

                for (int i = 0; i < nvalue; i++)
                {
                    for (int j = 0; j < nvalue; j++)
                    {
                        sum += m[i, j] * (x[i] * x[i] + x[i] + 1.0) * (x[j] * x[j] + x[j] + 1.0);
                    }
                }

                return sum;
            };

            double[][] b = new double[8][];
            b[0] = new double[16] { 0.22, 0.2, 0.19, 0.25, 0.15, 0.11, 0.12, 0.13, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            b[1] = new double[16] { -1.46, 0.0, -1.3, 1.82, -1.15, 0.0, 0.8, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            b[2] = new double[16] { 1.29, -0.89, 0.0, 0.0, -1.16, -0.96, 0.0, -0.49, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            b[3] = new double[16] { -1.1, -1.06, 0.95, -0.54, 0.0, -1.78, -0.41, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
            b[4] = new double[16] { 0.0, 0.0, 0.0, -1.43, 1.51, 0.59, -0.33, -0.43, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 };
            b[5] = new double[16] { 0.0, -1.72, -0.33, 0.0, 1.62, 1.24, 0.21, -0.26, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 };
            b[6] = new double[16] { 1.12, 0.0, 0.0, 0.31, 0.0, 0.0, 1.12, 0.0, -0.36, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 };
            b[7] = new double[16] { 0.0, 0.45, 0.26, -1.10, 0.58, 0.0, -1.03, 0.10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
            
            double[] c = new double[8] { 2.5, 1.1, -3.1, -3.5, 1.3, 2.1, 2.3, -1.5 };
            
            List<Func<double[], double>> eqConstraint = new List<Func<double[], double>>();

            for (int i = 0; i < nvalue/2; i++)
            {
                int k = i;
                Func<double[], double> eqc = (x) =>
                {
                    double sum = 0.0;
                    for (int j = 0; j < nvalue; j++)
                    {
                        sum += b[k][j] * x[j];
                    }
                    return sum - c[k];
                };

                eqConstraint.Add(eqc);
            }
            

            double[] upperBound = new double[nvalue];
            double[] lowerBound = new double[nvalue];
            double[] startValue = new double[nvalue];

            for (int i = 0; i < nvalue; i++)
            {
                lowerBound[i] = 0.0;
                upperBound[i] = 5.0;
                startValue[i] = 10.0;
            }
            double[] solution = new double[] { 0.03984735, 0.7919832, 0.2028703, 0.8443579, 1.126991, 0.9347387, 1.681962, 0.1553009, 1.56787, 0.0, 0.0, 0.0, 0.6602041, 0.0, 0.6742559, 0.0 };


            double fval = f(startValue);
            double fvalSol = f(solution);

            foreach (var item in eqConstraint)
            {
                double test = item(solution);
            }

            var res = quadraticProgramming.Minimize(f, eqConstraint, null,lowerBound, upperBound, startValue, 7000);

            double fval1 = f(res);


        }

        public static double[] StartSolution(
            List<Func<double[], double>> IneqConstraint,
            List<Func<double[], double>> EqConstraint,
            double[] startValue)
        {
            SQP quadraticProgramming = new SQP();

            Func<double[], double> f = (x) =>
            {
                return x[0];
            };
            
            List<Func<double[], double>> ineqVar = new List<Func<double[], double>>();

            foreach (var item in IneqConstraint)
            {
                Func<double[], double> buf = (x) =>
                {
                    return item(startValue) - x[0];
                };
                ineqVar.Add(buf);
            }

            var res =  quadraticProgramming.Minimize(f, ineqVar, null, new double[] { 0 }, 100);


            return res;
        }

        public static double[] GetInitialGuess(
            List<Func<double[], double>> IneqConstraint,
            List<Func<double[], double>> EqConstraint,
            double[] lowerBounds,
            double[] upperBounds,
            int count)
        {
            
            List<Func<double[], double>> allConstraint = new List<Func<double[], double>>(IneqConstraint);

            if(EqConstraint != null)
                allConstraint.AddRange(EqConstraint);

            SQP bfgs = new SQP();

            double[] startValue = new double[count];
            double[] fsolution = new double[count];
            double checkValue = double.MaxValue;
            int indexA = 0;
            int indexB = 0;

            ////Verifico che esista uno spigolo
            for (int i = 0; i < allConstraint.Count; i++)
            {
                for (int j = i + 1; j < allConstraint.Count; j++)
                {
                    Func<double[], double> f = (x) =>
                    {
                        double val = allConstraint[i](x) - allConstraint[j](x);
                        return val * val;
                    };

                    var solution = bfgs.Minimize(f, null, null, lowerBounds, upperBounds, startValue, 8);

                    //DEBUG Test solution
                    var t0 = f(solution);
                    //var t1 = allConstraint[i](solution);
                    //var t2 = allConstraint[j](solution);

                    if (t0 < checkValue)
                    {
                        bool check = true;
                        for (int k = 0; k < IneqConstraint.Count; k++)
                        {
                            if (IneqConstraint[k](solution) > 0)
                            {
                                check = false;
                                break;
                            }
                        }

                        if (check)
                        {
                            fsolution = solution;
                            checkValue = t0;
                            indexA = i;
                            indexB = j;
                        }
                    }
                }
            }

            //Reevaluate guess solution
            Func<double[], double> fun = (x) =>
            {
                double val = allConstraint[indexA](x) - allConstraint[indexB](x);
                return val * val;
            };

            fsolution = bfgs.Minimize(fun, null, null, lowerBounds, upperBounds, fsolution, 30);

            return fsolution;

        }

        private static List<Func<double[], double>> CreateBoundsConstraints(
            double[] lowerBound,
            double[] upperBound)
        {
            List<Func<double[], double>> boundsConstraints = new List<Func<double[], double>>();


            if (lowerBound != null)
            {
                for (int i = 0; i < lowerBound.Length; i++)
                {
                    int k = i;

                    Func<double[], double> boundConstraint = (x) =>
                    {
                        return -x[k] + lowerBound[k];
                    };

                    boundsConstraints.Add(boundConstraint);
                }
            }

            if (upperBound != null)
            {
                for (int i = 0; i < upperBound.Length; i++)
                {
                    int k = i;

                    Func<double[], double> boundConstraint = (x) =>
                    {
                        return x[k] - upperBound[k];
                    };

                    boundsConstraints.Add(boundConstraint);
                }
            }

            return boundsConstraints;
        }
    }
}
