using System;

namespace OptimizationMethods.Optimization
{
    public static class OptimizationHelper
    {
        private const double h = 10e-10;
        private const double h2 = h * 2;

        public static OptVector Derivative(Func<double[], double> f, OptVector x)
        {
            double[] result = new double[x.Count];

            for (int i = 0; i < x.Count; i++)
            {
                OptVector v = new OptVector(x);

                result[i] = (f(OptVector.Min(v, h2, i).MinArray) - 8.0 * f(OptVector.Min(v, h, i).MinArray) + 8.0 * f(OptVector.Sum(v, h, i).MinArray) - f(OptVector.Sum(v, h2, i).MinArray)) /
                            (h2 * 6.0);
            }

            return new OptVector(result);
        }

        public static OptVector Derivative(Func<double[], double>[] df, OptVector x)
        {
            double[] result = new double[x.Count];

            for (int i = 0; i < x.Count; i++)
            {
                OptVector v = new OptVector(x);

                result[i] = df[i](v.MinArray);
            }

            return new OptVector(result);
        }

        public static OptVector[] GetIdentity(int dim)
        {
            OptVector[] result = new OptVector[dim];

            for (int i = 0; i < dim; i++)
            {
                double[] v = new double[dim];
                v[i] = 1;
                result[i] = new OptVector(v);
            }

            return result;
        }
    }
}
