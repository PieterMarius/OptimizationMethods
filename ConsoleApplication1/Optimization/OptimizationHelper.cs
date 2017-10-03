using System;

namespace ConsoleApplication1.Optimization
{
    public static class OptimizationHelper
    {
        private const double h = 10e-10;
        private const double h2 = h * 2;

        public static MinVector Derivative(Func<double[], double> f, MinVector x)
        {
            double[] result = new double[x.Count];

            for (int i = 0; i < x.Count; i++)
            {
                MinVector v = new MinVector(x);

                result[i] = (f(MinVector.Min(v, h2, i).MinArray) - 8 * f(MinVector.Min(v, h, i).MinArray) + 8 * f(MinVector.Sum(v, h, i).MinArray) - f(MinVector.Sum(v, h2, i).MinArray)) /
                            (h2 * 6);
            }

            return new MinVector(result);
        }

        public static MinVector Derivative(Func<double[], double>[] df, MinVector x)
        {
            double[] result = new double[x.Count];

            for (int i = 0; i < x.Count; i++)
            {
                MinVector v = new MinVector(x);

                result[i] = df[i](v.MinArray);
            }

            return new MinVector(result);
        }

        public static MinVector[] GetIdentity(int dim)
        {
            MinVector[] result = new MinVector[dim];

            for (int i = 0; i < dim; i++)
            {
                double[] v = new double[dim];
                v[i] = 1;
                result[i] = new MinVector(v);
            }

            return result;
        }
    }
}
