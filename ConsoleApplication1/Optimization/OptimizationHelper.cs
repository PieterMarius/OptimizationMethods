using System;

namespace ConsoleApplication1.Optimization
{
    public static class OptimizationHelper
    {
        public const double h = 10e-10;

        public const double h2 = h * 2;

        public static double StrongWolfeLineSearch(
            Func<Vector, double> f,
            Vector d,
            Vector x0,
            double alpham)
        {
            int maxIter = 5;
            
            double alphap = 0;

            double c1 = 1E-4;
            double c2 = 0.5;

            Random rnd = new Random();

            double alphax = 1;

            double fx0 = f(x0);
            double gx0 = Derivative(f, x0) * d;

            double fxp = fx0;
            double gxp = gx0;

            for (int i = 1; i < maxIter; i++)
            {
                Vector xx = x0 + alphax * d;
                double fxx = f(xx);
                double gxx = Derivative(f, xx) * d;

                if (fxx > fx0 + c1 * alphax * gx0 ||
                    (i > 1 && fxx >= fxp))
                {
                    return Zoom(f, x0, d, alphap, alphax);
                }

                if (Math.Abs(gxx) <= -c2 * gx0)
                    return alphax;

                if (gxx >= 0)
                    return Zoom(f, x0, d, alphax, alphap);

                alphap = alphax;
                fxp = fxx;
                gxp = gxx;

                double r = 0.8; // rnd.NextDouble()

                alphax = alphax + (alpham - alphax) * r;
            }

            return alphax;
        }

        public static double Zoom(
            Func<Vector, double> f,
            Vector x0,
            Vector d,
            double alphal,
            double alphah)
        {
            int maxIter = 8;
            double c1 = 1E-4;
            double c2 = 0.5;

            double fx0 = f(x0);
            double gx0 = Derivative(f, x0) * d;

            double alphax = 0.0;

            for (int i = 0; i < maxIter; i++)
            {
                alphax = 0.5 * (alphal + alphah);
                Vector xx = x0 + alphax * d;
                double fxx = f(xx);
                double gxx = Derivative(f, xx) * d;
                Vector xl = x0 + alphal * d;
                double fxl = f(xl);

                if (fxx > fx0 + c1 * alphax * gx0 ||
                    fxx >= fxl)
                    alphah = alphax;
                else
                {
                    if (Math.Abs(gxx) <= -c2 * gx0)
                    {
                        return alphax;
                    }
                    if (gxx * (alphah - alphal) >= 0.0)
                        alphah = alphal;
                    alphal = alphax;
                }
            }
            return alphax;
        }

        public static double StrongWolfeLineSearch(
            Func<Vector, double> f,
            Func<Vector, double>[] df,
            Vector d,
            Vector x0,
            double alpham)
        {
            int maxIter = 5;

            double alphap = 0;

            double c1 = 1E-4;
            double c2 = 0.5;

            Random rnd = new Random();

            double alphax = 1;

            double fx0 = f(x0);
            double gx0 = Derivative(df, x0) * d;

            double fxp = fx0;
            double gxp = gx0;

            for (int i = 1; i < maxIter; i++)
            {
                Vector xx = x0 + alphax * d;
                double fxx = f(xx);
                double gxx = Derivative(df, xx) * d;

                if (fxx > fx0 + c1 * alphax * gx0 ||
                    (i > 1 && fxx >= fxp))
                {
                    return Zoom(f, df, x0, d, alphap, alphax);
                }

                if (Math.Abs(gxx) <= -c2 * gx0)
                    return alphax;

                if (gxx >= 0)
                    return Zoom(f, df, x0, d, alphax, alphap);

                alphap = alphax;
                fxp = fxx;
                gxp = gxx;

                alphax = alphax + (alpham - alphax) * rnd.NextDouble();
            }

            return alphax;
        }

        private static double Zoom(
            Func<Vector, double> f,
            Func<Vector, double>[] df,
            Vector x0,
            Vector d,
            double alphal,
            double alphah)
        {
            int maxIter = 8;
            double c1 = 1E-4;
            double c2 = 0.5;

            double fx0 = f(x0);
            double gx0 = Derivative(df, x0) * d;

            double alphax = 0.0;

            for (int i = 0; i < maxIter; i++)
            {
                alphax = 0.5 * (alphal + alphah);
                Vector xx = x0 + alphax * d;
                double fxx = f(xx);
                double gxx = Derivative(df, xx) * d;
                Vector xl = x0 + alphal * d;
                double fxl = f(xl);

                if (fxx > fx0 + c1 * alphax * gx0 ||
                    fxx >= fxl)
                    alphah = alphax;
                else
                {
                    if (Math.Abs(gxx) <= -c2 * gx0)
                    {
                        return alphax;
                    }
                    if (gxx * (alphah - alphal) >= 0.0)
                        alphah = alphal;
                    alphal = alphax;
                }
            }
            return alphax;
        }

        public static Vector Derivative(Func<Vector, double> f, Vector x)
        {
            double[] result = new double[x.Vars.Length];

            for (int i = 0; i < x.Vars.Length; i++)
            {
                Vector v = new Vector(x.Vars);

                result[i] = (f(Vector.Min(v, h2, i)) - 8 * f(Vector.Min(v, h, i)) + 8 * f(Vector.Sum(v, h, i)) - f(Vector.Sum(v, h2, i))) /
                            (h2 * 6);
            }

            return new Vector(result);
        }

        public static Vector Derivative(Func<Vector, double>[] df, Vector x)
        {
            double[] result = new double[x.Vars.Length];

            for (int i = 0; i < x.Vars.Length; i++)
            {
                Vector v = new Vector(x.Vars);

                result[i] = df[i](v);
            }

            return new Vector(result);
        }

        public static Vector[] GetIdentity(int dim)
        {
            Vector[] result = new Vector[dim];

            for (int i = 0; i < dim; i++)
            {
                double[] v = new double[dim];
                v[i] = 1;
                result[i] = new Vector(v);
            }

            return result;
        }
    }
}
