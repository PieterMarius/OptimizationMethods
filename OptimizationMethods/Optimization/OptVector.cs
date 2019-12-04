using System;
using System.Collections.Generic;
using System.Linq;

namespace OptimizationMethods.Optimization
{
    public struct OptVector
    {
        private double[] minArray { get; set; }

        public OptVector(double[] variables)
        {
            minArray = new double[variables.Length];
            Array.Copy(variables, minArray, variables.Length);
        }

        public OptVector(OptVector v)
        {
            minArray = new double[v.minArray.Length];
            Array.Copy(v.minArray, minArray, v.minArray.Length);
        }

        public OptVector(double var)
        {
            minArray = new double[1];
            minArray[0] = var;
        }

        public OptVector(int dim)
        {
            minArray = new double[dim];
        }

        public OptVector(int dim, double value)
        {
            minArray = new double[dim];

            for (int i = 0; i < minArray.Length; i++)
                minArray[i] = value;
        }

        public bool IsNull()
        {
            if (minArray != null)
                return false;

            return true;
        }

        public double this[int i]
        {
            get { return minArray[i]; }
            set { minArray[i] = value; }
        }

        public int Count
        {
            get { return minArray.Length; }
        }

        public double[] MinArray
        {
            get { return minArray; }
        }

        public static OptVector Populate(OptVector arr, double value)
        {
            for (int i = 0; i < arr.Count; i++)
            {
                arr.minArray[i] = value;
            }

            return arr;
        }

        public void Add(OptVector v)
        {
            var buf = minArray.ToList();

            buf.AddRange(v.minArray.ToList());

            minArray = buf.ToArray();
        }

        public void Add(double v)
        {
            var buf = minArray.ToList();

            buf.Add(v);

            minArray = buf.ToArray();
        }

        public void Add(List<double> v)
        {
            var buf = minArray.ToList();

            buf.AddRange(v);

            minArray = buf.ToArray();
        }

        public static OptVector Add(OptVector v, OptVector d)
        {
            var buf = v.minArray.ToList();

            buf.AddRange(d.minArray.ToList());

            return new OptVector(buf.ToArray());
        }

        public static OptVector operator +(OptVector a, OptVector b)
        {
            double[] result = new double[a.minArray.Length];
            for (int i = 0; i < a.minArray.Length; i++)
            {
                result[i] = a.minArray[i] + b.minArray[i];
            }

            return new OptVector(result);
        }

        public static OptVector operator +(OptVector a, double b)
        {
            double[] result = new double[a.minArray.Length];
            for (int i = 0; i < a.minArray.Length; i++)
            {
                result[i] = a.minArray[i] + b;
            }

            return new OptVector(result);
        }

        public static OptVector operator -(OptVector a, OptVector b)
        {
            double[] result = new double[a.minArray.Length];
            for (int i = 0; i < a.minArray.Length; i++)
            {
                result[i] = a.minArray[i] - b.minArray[i];
            }

            return new OptVector(result);
        }

        public static OptVector operator -(OptVector a, double b)
        {
            double[] result = new double[a.minArray.Length];
            for (int i = 0; i < a.minArray.Length; i++)
            {
                result[i] = a.minArray[i] - b;
            }

            return new OptVector(result);
        }

        public static OptVector operator *(OptVector a, double b)
        {
            double[] result = new double[a.minArray.Length];
            for (int i = 0; i < a.minArray.Length; i++)
            {
                result[i] = a.minArray[i] * b;
            }

            return new OptVector(result);
        }

        public static OptVector operator *(double b, OptVector a)
        {
            double[] result = new double[a.minArray.Length];
            for (int i = 0; i < a.minArray.Length; i++)
            {
                result[i] = a.minArray[i] * b;
            }

            return new OptVector(result);
        }


        public static double operator *(OptVector a, OptVector b)
        {
            double result = 0.0;
            for (int i = 0; i < a.minArray.Length; i++)
            {
                result += a.minArray[i] * b.minArray[i];
            }

            return result;
        }

        public static OptVector operator /(OptVector a, double b)
        {
            double[] result = new double[a.minArray.Length];
            for (int i = 0; i < a.minArray.Length; i++)
                result[i] = a.minArray[i] / b;

            return new OptVector(result);
        }

        public static OptVector operator /(double a, OptVector b)
        {
            double[] result = new double[b.minArray.Length];
            for (int i = 0; i < b.minArray.Length; i++)
                result[i] = a / b.minArray[i];

            return new OptVector(result);
        }

        public static OptVector[] Div(OptVector[] a, double b)
        {
            OptVector[] result = new OptVector[a.Length];

            for (int i = 0; i < a.Length; i++)
                result[i] = a[i] / b;

            return result;
        }

        public static OptVector[] Mult(OptVector a, OptVector b)
        {
            OptVector[] result = new OptVector[a.minArray.Length];

            for (int i = 0; i < a.minArray.Length; i++)
            {
                double[] v = new double[a.minArray.Length];

                for (int j = 0; j < a.minArray.Length; j++)
                    v[j] = a.minArray[i] * b.minArray[j];

                result[i] = new OptVector(v);
            }

            return result;
        }

        public static OptVector[] Mult(OptVector[] a, OptVector[] b)
        {
            OptVector[] result = new OptVector[a.Length];

            for (int i = 0; i < a.Length; i++)
            {
                double[] v = new double[a.Length];

                for (int j = 0; j < a.Length; j++)
                    v[j] = a[i] * GetColumn(b, j);

                result[i] = new OptVector(v);
            }

            return result;
        }

        public static OptVector Mult(OptVector[] a, OptVector b)
        {
            double[] result = new double[a.Length];

            for (int i = 0; i < a.Length; i++)
                result[i] = a[i] * b;

            return new OptVector(result);
        }

        public static OptVector Mult(OptVector a, OptVector[] b)
        {
            double[] result = new double[b.Length];

            for (int i = 0; i < b.Length; i++)
            {
                double sum = 0.0;
                for (int j = 0; j < b[i].Count; j++)
                {
                    sum += a[j] * b[j][i];
                }
                result[i] = sum;
            }

            return new OptVector(result);
        }

        public static OptVector[] SubtractFromIdentity(OptVector[] v)
        {
            OptVector[] result = new OptVector[v.Length];

            for (int i = 0; i < v.Length; i++)
            {
                double[] b = new double[v.Length];

                for (int j = 0; j < v.Length; j++)
                {
                    b[j] = (i == j) ? 1 - v[i].minArray[j] : -v[i].minArray[j];
                }

                result[i] = new OptVector(b);
            }

            return result;
        }

        public static OptVector[] InvertDiag(OptVector[] a)
        {
            OptVector[] result = new OptVector[a.Length];
            for (int i = 0; i < a.Length; i++)
            {
                result[i] = new OptVector(a.Length);
                if (a[i][i] != 0.0)
                    result[i][i] = 1.0 / a[i][i];
                else
                    result[i][i] = 0.0;
            }

            return result;
        }

        public static OptVector[] IncompleteCholesky(OptVector[] a)
        {
            OptVector[] result = new OptVector[a.Length];

            for (int i = 0; i < a.Length; i++)
            {
                result[i] = new OptVector(a.Length);
            }

            for (int i = 0; i < a.Length; i++)
            {
                double sum = 0.0;
                for (int k = 0; k < i - 1; k++)
                {
                    sum = a[i][k] * a[i][k];
                }

                result[i][i] = Math.Sqrt(a[i][i] - sum);

                for (int j = i + 1; j < a.Length; j++)
                {
                    sum = 0.0;
                    for (int k = 0; k < i - 1; k++)
                    {
                        sum = a[i][k] * a[j][k];
                    }
                    if(result[i][i] != 0.0)
                        result[j][i] = (1.0 / result[i][i]) * (a[j][i] - sum);
                }
            }

            return result;
        }

        public static OptVector Abs(OptVector a)
        {
            double[] result = new double[a.minArray.Length];
            for (int i = 0; i < a.minArray.Length; i++)
            {
                result[i] = Math.Abs(a.minArray[i]);
            }

            return new OptVector(result);
        }

        public static OptVector Min(OptVector a, double v, int i)
        {
            a.minArray[i] = a.minArray[i] - v;

            return new OptVector(a.minArray);
        }

        public static OptVector Sum(OptVector a, double v, int i)
        {
            a.minArray[i] = a.minArray[i] + v;

            return new OptVector(a.minArray);
        }

        public static OptVector[] Sum(OptVector[] a, OptVector[] b)
        {
            OptVector[] result = new OptVector[a.Length];
            for (int i = 0; i < a.Length; i++)
            {
                result[i] = a[i] + b[i];
            }

            return result;
        }

        public static OptVector[] Sum(OptVector[] a, double b)
        {
            OptVector[] result = new OptVector[a.Length];
            for (int i = 0; i < a.Length; i++)
            {
                result[i] = a[i] + b;
            }

            return result;
        }

        public static OptVector GetColumn(OptVector[] a, int index)
        {
            double[] result = new double[a.Length];

            for (int i = 0; i < a.Length; i++)
                result[i] = a[i].minArray[index];

            return new OptVector(result);
        }

        public static OptVector GetIndex(OptVector a, int index)
        {
            double[] vt = new double[a.minArray.Length];

            vt[index] = a.minArray[index];

            return new OptVector(vt);
        }

        public static OptVector[] Transpose(OptVector[] src)
        {
            OptVector[] result = new OptVector[src.Length];

            for (int i = 0; i < result.Length; i++)
                result[i] = new OptVector(result.Length);

            for (int i = 0; i < src.Length; i++)
            {
                OptVector buf = new OptVector(src.Length);

                for (int j = 0; j < src.Length; j++)
                    buf.minArray[j] = src[j].minArray[i];

                result[i] = buf;

            }

            return result;
        }

        public static bool CheckPositiveDefiniteMatrix(OptVector[] m)
        {
            for (int i = 0; i < m.Length; i++)
            {
                double sum = 0.0;

                if (m[i].minArray[i] < 0.0)
                    return false;

                for (int j = 0; j < m[i].minArray.Length; j++)
                    sum += Math.Abs(m[i].minArray[j]);

                if (m[i].minArray[i] < sum - m[i].minArray[i])
                    return false;

            }

            return true;

        }

        public static bool Equals(OptVector a, OptVector b)
        {
            if (a.Count != b.Count)
                return false;

            for (int i = 0; i < a.Count; i++)
            {
                if (a[i] != b[i])
                    return false;
            }

            return true;
        }

        public static bool Equals(OptVector[] a, OptVector[] b)
        {
            if (a.Length != b.Length)
                return false;

            for (int i = 0; i < a.Length; i++)
            {
                if (!Equals(a[i], b[i]))
                    return false;
            }

            return true;
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

        public static OptVector[] GetIdentity(int dim, OptVector values)
        {
            OptVector[] result = new OptVector[dim];

            for (int i = 0; i < dim; i++)
            {
                double[] v = new double[dim];
                v[i] = values[i];
                result[i] = new OptVector(v);
            }

            return result;
        }

        public double Length()
        {
            return Math.Sqrt(this * this);
        }

        public OptVector Normalize()
        {
            return this / this.Length();
        }
    }
}
