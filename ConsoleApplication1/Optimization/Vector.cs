using System;
using System.Collections.Generic;
using System.Linq;

namespace ConsoleApplication1.Optimization
{
    public struct Vector
    {
        public double[] Vars { get; set; }
        
        public Vector(double[] variables)
        {
            Vars = new double[variables.Length];
            Array.Copy(variables, Vars, variables.Length);
        }

        public Vector(double var)
        {
            Vars = new double[1];
            Vars[0] = var;
        }

        public Vector(int dim)
        {
            Vars = new double[dim];
        }

        public double this[int i]
        {
            get { return Vars[i]; }
            set { Vars[i] = value; }
        }

        public int Count()
        {
            return Vars.Length;
        }

        public static Vector Populate(Vector arr, double value)
        {
            for (int i = 0; i < arr.Count(); i++)
            {
                arr.Vars[i] = value;
            }

            return arr;
        }

        public void Add(Vector v)
        {
            var buf = Vars.ToList();

            buf.AddRange(v.Vars.ToList());

            Vars = buf.ToArray();
        }

        public void Add(double v)
        {
            var buf = Vars.ToList();

            buf.Add(v);

            Vars = buf.ToArray();
        }

        public void Add(List<double> v)
        {
            var buf = Vars.ToList();

            buf.AddRange(v);

            Vars = buf.ToArray();
        }

        public static Vector Add(Vector v, Vector d)
        {
            var buf = v.Vars.ToList();

            buf.AddRange(d.Vars.ToList());

            return new Vector(buf.ToArray());
        }

        public static Vector operator +(Vector a, Vector b)
        {
            double[] result = new double[a.Vars.Length];
            for (int i = 0; i < a.Vars.Length; i++)
            {
                result[i] = a.Vars[i] + b.Vars[i];
            }

            return new Vector(result);
        }

        public static Vector operator +(Vector a, double b)
        {
            double[] result = new double[a.Vars.Length];
            for (int i = 0; i < a.Vars.Length; i++)
            {
                result[i] = a.Vars[i] + b;
            }

            return new Vector(result);
        }

        public static Vector operator -(Vector a, Vector b)
        {
            double[] result = new double[a.Vars.Length];
            for (int i = 0; i < a.Vars.Length; i++)
            {
                result[i] = a.Vars[i] - b.Vars[i];
            }

            return new Vector(result);
        }

        public static Vector operator -(Vector a, double b)
        {
            double[] result = new double[a.Vars.Length];
            for (int i = 0; i < a.Vars.Length; i++)
            {
                result[i] = a.Vars[i] - b;
            }

            return new Vector(result);
        }

        public static Vector operator *(Vector a, double b)
        {
            double[] result = new double[a.Vars.Length];
            for (int i = 0; i < a.Vars.Length; i++)
            {
                result[i] = a.Vars[i] * b;
            }

            return new Vector(result);
        }

        public static Vector operator *(double b, Vector a)
        {
            double[] result = new double[a.Vars.Length];
            for (int i = 0; i < a.Vars.Length; i++)
            {
                result[i] = a.Vars[i] * b;
            }

            return new Vector(result);
        }

        
        public static double operator *(Vector a, Vector b)
        {
            double result = 0.0;
            for (int i = 0; i < a.Vars.Length; i++)
            {
                result += a.Vars[i] * b.Vars[i];
            }

            return result;
        }

        public static Vector operator /(Vector a, double b)
        {
            double[] result = new double[a.Vars.Length];
            for (int i = 0; i < a.Vars.Length; i++)
                result[i] = a.Vars[i] / b;

            return new Vector(result);
        }

        public static Vector operator /(double a, Vector b)
        {
            double[] result = new double[b.Vars.Length];
            for (int i = 0; i < b.Vars.Length; i++)
                result[i] = a / b.Vars[i];

            return new Vector(result);
        }

        public static Vector[] Div(Vector[] a, double b)
        {
            Vector[] result = new Vector[a.Length];

            for (int i = 0; i < a.Length; i++)
                result[i] = a[i] / b;
            
            return result;
        }

        public static Vector[] Mult(Vector a, Vector b)
        {
            Vector[] result = new Vector[a.Vars.Length];

            for (int i = 0; i < a.Vars.Length; i++)
            {
                double[] v = new double[a.Vars.Length];

                for (int j = 0; j < a.Vars.Length; j++)
                    v[j] = a.Vars[i] * b.Vars[j];
                
                result[i] = new Vector(v);
            }

            return result;
        }

        public static Vector[] Mult(Vector[] a, Vector[] b)
        {
            Vector[] result = new Vector[a.Length];

            for (int i = 0; i < a.Length; i++)
            {
                double[] v = new double[a.Length];

                for (int j = 0; j < a.Length; j++)
                    v[j] = a[i] * GetColumn(b, j);
                                
                result[i] = new Vector(v);
            }

            return result;
        }

        public static Vector Mult(Vector[] a, Vector b)
        {
            double[] result = new double[a.Length];

            for (int i = 0; i < a.Length; i++)
                result[i] = a[i] * b;

            return new Vector(result);
        }

        public static Vector Mult(Vector a, Vector[] b)
        {
            double[] result = new double[b.Length];

            for (int i = 0; i < b.Length; i++)
            {
                double sum = 0.0;
                for (int j = 0; j < b[i].Count(); j++)
                {
                    sum += a[j] * b[j][i];
                }
                result[i] = sum;
            }

            return new Vector(result);
        }
                       
        public static Vector[] SubtractFromIdentity(Vector[] v)
        {
            Vector[] result = new Vector[v.Length];

            for (int i = 0; i < v.Length; i++)
            {
                double[] b = new double[v.Length];

                for (int j = 0; j < v.Length; j++)
                {
                    if (i == j)
                        b[j] = 1 - v[i].Vars[j];
                    else
                        b[j] = -v[i].Vars[j];
                }

                result[i] = new Vector(b);
            }

            return result;
        }
                
        public static Vector Abs(Vector a)
        {
            double[] result = new double[a.Vars.Length];
            for (int i = 0; i < a.Vars.Length; i++)
            {
                result[i] = Math.Abs(a.Vars[i]);
            }

            return new Vector(result);
        }

        public static Vector Min(Vector a, double v, int i)
        {
            a.Vars[i] = a.Vars[i] - v;

            return new Vector(a.Vars);
        }

        public static Vector Sum(Vector a, double v, int i)
        {
            a.Vars[i] = a.Vars[i] + v;

            return new Vector(a.Vars);
        }

        public static Vector[] Sum(Vector[] a, Vector[] b)
        {
            Vector[] result = new Vector[a.Length];
            for (int i = 0; i < a.Length; i++)
            {
                result[i] = a[i] + b[i];
            }

            return result;
        }

        public static Vector GetColumn(Vector[] a, int index)
        {
            double[] result = new double[a.Length];

            for (int i = 0; i < a.Length; i++)
                result[i] = a[i].Vars[index];

            return new Vector(result);
        }

        public static Vector GetIndex(Vector a, int index)
        {
            double[] vt = new double[a.Vars.Length];

            vt[index] = a.Vars[index];

            return new Vector(vt);
        }

        public static Vector[] Transpose(Vector[] src)
        {
            Vector[] result = new Vector[src.Length];
            for (int i = 0; i < result.Length; i++)
                result[i] = new Vector(result.Length);
            
            for(int i = 0; i< src.Length; i++)
            {
                Vector buf = new Vector(src.Length);
                for (int j = 0; j < src.Length; j++)
                    buf.Vars[j] = src[j].Vars[i];

                result[i] = buf;
            }

            return result;
        }

        public static bool Equals(Vector[] a, Vector[] b)
        {
            if (a.Length != b.Length)
                return false;

            for (int i = 0; i < a.Length; i++)
            {
                if (a[i].Vars.Length != b[i].Vars.Length)
                    return false;

                for (int j = 0; j < a[i].Vars.Length; j++)
                {
                    if (Math.Abs(a[i].Vars[j] - b[i].Vars[j]) > 1E-10)
                        return false;
                }
            }

            return true;
        }

        public static bool CheckPositiveMatrix(Vector[] m)
        {
            for (int i = 0; i < m.Length; i++)
            {
                if (m[i].Vars[i] < 0.0)
                    return false;
            }

            return true;
        }

        public static bool CheckPositiveDefiniteMatrix(Vector[] m)
        {
            for (int i = 0; i < m.Length; i++)
            {
                double sum = 0.0;

                if (m[i].Vars[i] < 0.0)
                    return false;

                for (int j = 0; j < m[i].Vars.Length; j++)
                {
                    sum += Math.Abs(m[i].Vars[j]);
                }

                if (m[i].Vars[i] < sum - m[i].Vars[i])
                    return false;
            }

            return true;
        }

        public static Vector[] SetIdentity(Vector a)
        {
            Vector[] result = OptimizationHelper.GetIdentity(a.Vars.Length);

            for (int i = 0; i < a.Vars.Length; i++)
                result[i].Vars[i] = a.Vars[i];

            return result;
        }
    }
}
