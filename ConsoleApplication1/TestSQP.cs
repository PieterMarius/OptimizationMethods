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

            Func<Vector, double> f = (x) =>
            {
                return Math.Pow(x.Vars[0], 2) +
                       Math.Pow(x.Vars[1], 2);
            };

            List<Func<Vector, double>> eqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> eqConstraint1 = (x) =>
            {
                return -x.Vars[0] -
                       x.Vars[1] + 2;
            };

            eqConstraint.Add(eqConstraint1);

            var res = quadraticProgramming.Minimize(f, eqConstraint, null, new double[] { 0, 0 }, 50);
        }

        public static void Test1()
        {
            SQP quadraticProgramming = new SQP();

            Func<Vector, double> f = (x) =>
            {
                return Math.Pow(x.Vars[0] - 1, 2) +
                       Math.Pow(x.Vars[1] - 3, 2);
            };

            List<Func<Vector, double>> eqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> eqConstraint1 = (x) =>
            {
                return -x.Vars[0] +
                       Math.Pow(x.Vars[1], 2) - 1;
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

            Func<Vector, double> f = (x) =>
            {
                return -(5 - Math.Pow(x.Vars[0] - 2, 2) -
                       2 * Math.Pow(x.Vars[1] - 1, 2));
            };

            List<Func<Vector, double>> eqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> eqConstraint1 = (x) =>
            {
                return -x.Vars[0] -
                       4 * x.Vars[1] + 3;
            };

            eqConstraint.Add(eqConstraint1);

            var res = quadraticProgramming.Minimize(f, eqConstraint, null, new double[] { 0, 0 }, 50);
        }

        /// <summary>
        /// x1= 1.0, x2 = 0.0, x3 = 0.0
        /// </summary>
        public static void Test3()
        {
            SQP quadraticProgramming = new SQP();

            Func<Vector, double> f = (x) =>
            {
                return x.Vars[0] + x.Vars[1] + Math.Pow(x.Vars[2], 2);
            };

            List<Func<Vector, double>> eqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> eqConstraint1 = (x) =>
            {
                return x.Vars[0] - 1;
            };

            Func<Vector, double> eqConstraint2 = (x) =>
            {
                return Math.Pow(x.Vars[0], 2) + Math.Pow(x.Vars[1], 2) - 1;
            };

            eqConstraint.Add(eqConstraint1);
            eqConstraint.Add(eqConstraint2);

            var res = quadraticProgramming.Minimize(f, eqConstraint, null, new double[] { 0, 0, 0 }, 1000);
        }

        public static void Test4()
        {
            SQP quadraticProgramming = new SQP();

            Func<Vector, double> f = (x) =>
            {

                return 400 * Math.Pow(x.Vars[0], 2) +
                       800 * Math.Pow(x.Vars[1], 2) +
                       200 * x.Vars[0] * x.Vars[1] +
                       1600 * Math.Pow(x.Vars[2], 2) +
                       400 * x.Vars[1] * x.Vars[2];
            };

            List<Func<Vector, double>> eqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> eqConstraint1 = (x) =>
            {
                return 1.2 - x.Vars[0] - x.Vars[1] - 1.5 * x.Vars[2];
            };

            Func<Vector, double> eqConstraint2 = (x) =>
            {
                return 1.0 - x.Vars[0] - x.Vars[1] - x.Vars[2];
            };

            eqConstraint.Add(eqConstraint2);
            eqConstraint.Add(eqConstraint1);

            var res = quadraticProgramming.Minimize(f, eqConstraint, null, new double[] { 0, 0, 0 }, 50);
        }

        public static void Test5()
        {
            SQP quadraticProgramming = new SQP();

            Func<Vector, double> f = (x) =>
            {

                return Math.Pow(x[0] - 2, 2) +
                       2 * Math.Pow(x[1] - 1, 2);
            };

            List<Func<Vector, double>> inqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> inqConstraint1 = (x) =>
            {
                return 3 - x[0] - 4 * x[1];
            };

            Func<Vector, double> inqConstraint2 = (x) =>
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

            Func<Vector, double> f = (x) =>
            {

                return -(x[0] * (30 - x[0]) + x[1] * (50 - 2 * x[1]) - 3 * x[0] - 5 * x[1] - 10 * x[2]);
            };

            List<Func<Vector, double>> inqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> inqConstraint1 = (x) =>
            {
                return x[0] + x[1] - x[2];
            };

            Func<Vector, double> inqConstraint2 = (x) =>
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

            Func<Vector, double> f = (x) =>
            {

                return -(Math.Sin(x[0]) * Math.Cos(x[1]) + Math.Cos(x[0]) * Math.Sin(x[1]));
            };

            List<Func<Vector, double>> inqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> inqConstraint1 = (x) =>
            {
                return -x[0] - x[1];
            };

            Func<Vector, double> inqConstraint2 = (x) =>
            {
                return x[0] + x[1] - Math.PI;
            };

            inqConstraint.Add(inqConstraint1);
            inqConstraint.Add(inqConstraint2);

            List<Func<Vector, double>> eqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> eqConstraint1 = (x) =>
            {
                return x[0] - Math.Pow(x[1], 3);
            };

            eqConstraint.Add(eqConstraint1);

            var res = quadraticProgramming.Minimize(f, eqConstraint, inqConstraint, new double[] { 0, 0 }, 100);
        }

        /// <summary>
        /// x1 = -1.0, x2 = -1.0
        /// </summary>
        public static void Test8()
        {
            SQP quadraticProgramming = new SQP();

            Func<Vector, double> f = (x) =>
            {

                return x[0] + x[1];
            };

            List<Func<Vector, double>> inqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> inqConstraint1 = (x) =>
            {
                return -2 + Math.Pow(x[0], 2) + Math.Pow(x[1], 2);
            };

            inqConstraint.Add(inqConstraint1);

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { 0, 0 }, 100);
        }


        /// <summary>
        /// Rosenbrock's function (x1 = 0.4149, x2 = 0.1701)
        /// </summary>
        public static void Test9()
        {
            SQP quadraticProgramming = new SQP();

            Func<Vector, double> f = (x) =>
            {

                return Math.Pow(1 - x[0], 2) + 100 * Math.Pow(x[1] - Math.Pow(x[0], 2), 2);
            };

            List<Func<Vector, double>> inqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> inqConstraint1 = (x) =>
            {
                return x[0] + 2 * x[1] - 1;
            };

            inqConstraint.Add(inqConstraint1);

            List<Func<Vector, double>> eqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> eqConstraint1 = (x) =>
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

            Func<Vector, double> f = (x) =>
            {

                return 1 + x[0] / (1 + x[1]) - 3 * x[0] * x[1] + x[1] * (1 + x[0]);
            };

            List<Func<Vector, double>> inqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> inqConstraint1 = (x) =>
            {
                return x[0] - 1;
            };

            Func<Vector, double> inqConstraint2 = (x) =>
            {
                return x[1] - 2;
            };

            Func<Vector, double> inqConstraint3 = (x) =>
            {
                return -x[0];
            };


            Func<Vector, double> inqConstraint4 = (x) =>
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

            Func<Vector, double> f = (x) =>
            {

                return Math.Pow(1.0 - x[0], 2) + 100.0 * Math.Pow(x[1] - Math.Pow(x[0], 2), 2);
            };

            List<Func<Vector, double>> inqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> eqConstraint1 = (x) =>
            {
                return Math.Pow(x[0] - (1.0 / 3.0), 2) + Math.Pow(x[1] - (1.0 / 3.0), 2) - Math.Pow(1.0 / 3.0, 2);
            };

            Func<Vector, double> inqConstraint1 = (x) =>
            {
                return x[0] - 0.5;
            };

            Func<Vector, double> inqConstraint2 = (x) =>
            {
                return x[1] - 0.8;
            };

            Func<Vector, double> inqConstraint3 = (x) =>
            {
                return -x[0];
            };

            Func<Vector, double> inqConstraint4 = (x) =>
            {
                return -x[1] + 0.2;
            };

            inqConstraint.Add(eqConstraint1);
            inqConstraint.Add(inqConstraint1);
            inqConstraint.Add(inqConstraint2);
            inqConstraint.Add(inqConstraint3);
            inqConstraint.Add(inqConstraint4);

            Vector upperBound = new Vector(new double[] { 0.5, 0.8 });
            Vector lowerBound = new Vector(new double[] { 0.0, 0.2 });

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, upperBound, lowerBound, new double[] { 0.25, 0.25 }, 1000);
        }


        /// <summary>
        /// x1 = 1.4, x2 = 1.7
        /// </summary>
        public static void Test12()
        {
            SQP quadraticProgramming = new SQP();

            Func<Vector, double> f = (x) =>
            {

                return Math.Pow(x[0] - 1.0, 2) + Math.Pow(x[1] - 2.5, 2);
            };

            List<Func<Vector, double>> inqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> inqConstraint1 = (x) =>
            {
                return -x[0] + 2 * x[1] - 2.0;
            };

            Func<Vector, double> inqConstraint2 = (x) =>
            {
                return x[0] + x[1] - 6.0;
            };

            Func<Vector, double> inqConstraint3 = (x) =>
            {
                return x[0] - 2.0 * x[1] - 2.0;
            };

            Func<Vector, double> inqConstraint4 = (x) =>
            {
                return -x[0];
            };

            Func<Vector, double> inqConstraint5 = (x) =>
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

            Func<Vector, double> f = (x) =>
            {

                return Math.Pow(1.0 - x[0], 2) + 100.0 * Math.Pow(x[1] - Math.Pow(x[0], 2), 2);
            };

            List<Func<Vector, double>> inqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> inqConstraint1 = (x) =>
            {
                return x[0] - 0.5;
            };

            Func<Vector, double> inqConstraint2 = (x) =>
            {
                return x[1] - 0.8;
            };

            Func<Vector, double> inqConstraint3 = (x) =>
            {
                return -x[0];
            };

            Func<Vector, double> inqConstraint4 = (x) =>
            {
                return -x[1] + 0.2;
            };

            Func<Vector, double> eqConstraint1 = (x) =>
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

            Func<Vector, double> f = (x) =>
            {

                return Math.Pow(1.0 - x[0], 2) + 100.0 * Math.Pow(x[1] - Math.Pow(x[0], 2), 2);
            };

            List<Func<Vector, double>> inqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> eqConstraint1 = (x) =>
            {
                return Math.Pow(x[0], 2) + Math.Pow(x[1], 2) - 1;
            };

            inqConstraint.Add(eqConstraint1);

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { 0.0, 0.0 }, 100);
        }

        public static void Test15()
        {
            SQP quadraticProgramming = new SQP();

            Func<Vector, double> f = (x) =>
            {

                return -(2 * x[0] * x[1] + 2 * x[0] - Math.Pow(x[0], 2) - 2 * Math.Pow(x[1], 2));
            };

            List<Func<Vector, double>> inqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> inqConstraint1 = (x) =>
            {
                return -x[1] + 1.0;
            };

            inqConstraint.Add(inqConstraint1);

            List<Func<Vector, double>> eqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> eqConstraint1 = (x) =>
            {
                return Math.Pow(x[0], 3) - x[1];
            };

            eqConstraint.Add(eqConstraint1);

            var res = quadraticProgramming.Minimize(f, eqConstraint, inqConstraint, new double[] { 0.0, 0.0 }, 100);
        }


        /// <summary>
        /// x1 = 2.0, x2 = 0.0
        /// </summary>
        public static void Test16()
        {
            SQP quadraticProgramming = new SQP();

            Func<Vector, double> f = (x) =>
            {
                return 0.01 * Math.Pow(x[0], 2) + Math.Pow(x[1], 2) - 100.0;
            };

            List<Func<Vector, double>> inqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> inqConstraint1 = (x) =>
            {
                return x[0] - 50;
            };

            Func<Vector, double> inqConstraint2 = (x) =>
            {
                return x[1] - 50;
            };

            Func<Vector, double> inqConstraint3 = (x) =>
            {
                return -x[0] + 2.0;
            };

            Func<Vector, double> inqConstraint4 = (x) =>
            {
                return -x[1] - 50;
            };

            Func<Vector, double> inqConstraint5 = (x) =>
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

            Func<Vector, double> f = (x) =>
            {

                return Math.Pow(1.0 - x[0], 2) + 100.0 * Math.Pow(x[1] - Math.Pow(x[0], 2), 2);
            };

            List<Func<Vector, double>> inqConstraint = new List<Func<Vector, double>>();

            Func<Vector, double> eqConstraint1 = (x) =>
            {
                return Math.Pow(x[0], 2) + Math.Pow(x[1], 2) - 1.5;
            };

            inqConstraint.Add(eqConstraint1);

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { 1.0, -0.5 }, 100);
        }

        public static void Test18()
        {
            SQP quadraticProgramming = new SQP();

            Func<Vector, double> f = (x) =>
            {

                return Math.Pow(x[0], 2) + Math.Pow(x[1], 2) + Math.Pow(x[2], 2);
            };
            //−x1 − x2 exp (x3y) − exp(2y) + 2 exp(4y) ≥ 0 for all y ∈ [0, 1]

            List <Func<Vector, double>> inqConstraint = new List<Func<Vector, double>>();


            int nConstraints = 100;
            double step = 1.0 / nConstraints;
            double aStep = 0.0;

            for (int i = 0; i < nConstraints; i++)
            {
                aStep += step;

                double inner = -aStep;

                Func<Vector, double> eqConstraint1 = (x) =>
                {
                    return x[0] + x[1] * Math.Exp(x[2] * inner) + Math.Exp(2 * inner) - 2 * Math.Exp(4 * inner);
                };

                inqConstraint.Add(eqConstraint1);
            }
                      

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { 1.0, -1.0, 2.0}, 10000);
        }

        /// <summary>
        /// x1 = 1.6667, x2 = 2.11111
        /// </summary>
        public static void Test19()
        {
            SQP quadraticProgramming = new SQP();

            Func<Vector, double> f = (x) =>
            {

                return Math.Pow(x[0], 4) + Math.Pow(x[1], 4);
            };
            

            List<Func<Vector, double>> inqConstraint = new List<Func<Vector, double>>();

                   

            Func<Vector, double> eqConstraint1 = (x) =>
            {
                return Math.Pow(x[0],2) - x[0] - x[1] + 1;
            };

            Func<Vector, double> eqConstraint2 = (x) =>
            {
                return Math.Pow(x[0], 2) - 4 * x[0] - x[1] + 6;
            };

            Func<Vector, double> eqConstraint3 = (x) =>
            {
                return Math.Pow(x[0], 2) - 3 * x[0] + x[1] - 2;
            };

            inqConstraint.Add(eqConstraint1);
            inqConstraint.Add(eqConstraint2);
            inqConstraint.Add(eqConstraint3);

            var res = quadraticProgramming.Minimize(f, null, inqConstraint, new double[] { - 4.0, -4.0 }, 10000);
        }


    }
}
