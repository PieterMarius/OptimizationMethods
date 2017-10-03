using System;
using System.Collections.Generic;
using System.Linq;

namespace ConsoleApplication1.Optimization
{
    public struct MinVector
    {
        private double[] minArray { get; set; }

        public MinVector(double[] variables)
        {
            minArray = new double[variables.Length];
            Array.Copy(variables, minArray, variables.Length);
        }

        public MinVector(MinVector v)
        {
            minArray = new double[v.minArray.Length];
            Array.Copy(v.minArray, minArray, v.minArray.Length);
        }

        public MinVector(double var)
        {
            minArray = new double[1];
            minArray[0] = var;
        }

        public MinVector(int dim)
        {
            minArray = new double[dim];
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

        public static MinVector Populate(MinVector arr, double value)
        {
            for (int i = 0; i < arr.Count; i++)
            {
                arr.minArray[i] = value;
            }

            return arr;
        }

        public void Add(MinVector v)
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

        public static MinVector Add(MinVector v, MinVector d)
        {
            var buf = v.minArray.ToList();

            buf.AddRange(d.minArray.ToList());

            return new MinVector(buf.ToArray());
        }

        public static MinVector operator +(MinVector a, MinVector b)
        {
            double[] result = new double[a.minArray.Length];
            for (int i = 0; i < a.minArray.Length; i++)
            {
                result[i] = a.minArray[i] + b.minArray[i];
            }

            return new MinVector(result);
        }

        public static MinVector operator +(MinVector a, double b)
        {
            double[] result = new double[a.minArray.Length];
            for (int i = 0; i < a.minArray.Length; i++)
            {
                result[i] = a.minArray[i] + b;
            }

            return new MinVector(result);
        }

        public static MinVector operator -(MinVector a, MinVector b)
        {
            double[] result = new double[a.minArray.Length];
            for (int i = 0; i < a.minArray.Length; i++)
            {
                result[i] = a.minArray[i] - b.minArray[i];
            }

            return new MinVector(result);
        }

        public static MinVector operator -(MinVector a, double b)
        {
            double[] result = new double[a.minArray.Length];
            for (int i = 0; i < a.minArray.Length; i++)
            {
                result[i] = a.minArray[i] - b;
            }

            return new MinVector(result);
        }

        public static MinVector operator *(MinVector a, double b)
        {
            double[] result = new double[a.minArray.Length];
            for (int i = 0; i < a.minArray.Length; i++)
            {
                result[i] = a.minArray[i] * b;
            }

            return new MinVector(result);
        }

        public static MinVector operator *(double b, MinVector a)
        {
            double[] result = new double[a.minArray.Length];
            for (int i = 0; i < a.minArray.Length; i++)
            {
                result[i] = a.minArray[i] * b;
            }

            return new MinVector(result);
        }


        public static double operator *(MinVector a, MinVector b)
        {
            double result = 0.0;
            for (int i = 0; i < a.minArray.Length; i++)
            {
                result += a.minArray[i] * b.minArray[i];
            }

            return result;
        }

        public static MinVector operator /(MinVector a, double b)
        {
            double[] result = new double[a.minArray.Length];
            for (int i = 0; i < a.minArray.Length; i++)
                result[i] = a.minArray[i] / b;

            return new MinVector(result);
        }

        public static MinVector operator /(double a, MinVector b)
        {
            double[] result = new double[b.minArray.Length];
            for (int i = 0; i < b.minArray.Length; i++)
                result[i] = a / b.minArray[i];

            return new MinVector(result);
        }

        public static MinVector[] Div(MinVector[] a, double b)
        {
            MinVector[] result = new MinVector[a.Length];

            for (int i = 0; i < a.Length; i++)
                result[i] = a[i] / b;

            return result;
        }

        public static MinVector[] Mult(MinVector a, MinVector b)
        {
            MinVector[] result = new MinVector[a.minArray.Length];

            for (int i = 0; i < a.minArray.Length; i++)
            {
                double[] v = new double[a.minArray.Length];

                for (int j = 0; j < a.minArray.Length; j++)
                    v[j] = a.minArray[i] * b.minArray[j];

                result[i] = new MinVector(v);
            }

            return result;
        }

        public static MinVector[] Mult(MinVector[] a, MinVector[] b)
        {
            MinVector[] result = new MinVector[a.Length];

            for (int i = 0; i < a.Length; i++)
            {
                double[] v = new double[a.Length];

                for (int j = 0; j < a.Length; j++)
                    v[j] = a[i] * GetColumn(b, j);

                result[i] = new MinVector(v);
            }

            return result;
        }

        public static MinVector Mult(MinVector[] a, MinVector b)
        {
            double[] result = new double[a.Length];

            for (int i = 0; i < a.Length; i++)
                result[i] = a[i] * b;

            return new MinVector(result);
        }

        public static MinVector Mult(MinVector a, MinVector[] b)
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

            return new MinVector(result);
        }

        public static MinVector[] SubtractFromIdentity(MinVector[] v)
        {
            MinVector[] result = new MinVector[v.Length];

            for (int i = 0; i < v.Length; i++)
            {
                double[] b = new double[v.Length];

                for (int j = 0; j < v.Length; j++)
                {
                    b[j] = (i == j) ? 1 - v[i].minArray[j] : -v[i].minArray[j];
                }

                result[i] = new MinVector(b);
            }

            return result;
        }

        public static MinVector Abs(MinVector a)
        {
            double[] result = new double[a.minArray.Length];
            for (int i = 0; i < a.minArray.Length; i++)
            {
                result[i] = Math.Abs(a.minArray[i]);
            }

            return new MinVector(result);
        }

        public static MinVector Min(MinVector a, double v, int i)
        {
            a.minArray[i] = a.minArray[i] - v;

            return new MinVector(a.minArray);
        }

        public static MinVector Sum(MinVector a, double v, int i)
        {
            a.minArray[i] = a.minArray[i] + v;

            return new MinVector(a.minArray);
        }

        public static MinVector[] Sum(MinVector[] a, MinVector[] b)
        {
            MinVector[] result = new MinVector[a.Length];
            for (int i = 0; i < a.Length; i++)
            {
                result[i] = a[i] + b[i];
            }

            return result;
        }

        public static MinVector GetColumn(MinVector[] a, int index)
        {
            double[] result = new double[a.Length];

            for (int i = 0; i < a.Length; i++)
                result[i] = a[i].minArray[index];

            return new MinVector(result);
        }

        public static MinVector GetIndex(MinVector a, int index)
        {
            double[] vt = new double[a.minArray.Length];

            vt[index] = a.minArray[index];

            return new MinVector(vt);
        }

        public static MinVector[] Transpose(MinVector[] src)
        {
            MinVector[] result = new MinVector[src.Length];

            for (int i = 0; i < result.Length; i++)
                result[i] = new MinVector(result.Length);

            for (int i = 0; i < src.Length; i++)
            {
                MinVector buf = new MinVector(src.Length);

                for (int j = 0; j < src.Length; j++)
                    buf.minArray[j] = src[j].minArray[i];

                result[i] = buf;

            }

            return result;
        }

        public static bool CheckPositiveDefiniteMatrix(MinVector[] m)
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

        public static bool Equals(MinVector a, MinVector b)
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

        public static bool Equals(MinVector[] a, MinVector[] b)
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

        public double Length()
        {
            return Math.Sqrt(this * this);
        }
    }
}
