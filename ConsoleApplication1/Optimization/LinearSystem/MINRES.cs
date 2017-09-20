using System;

namespace ConsoleApplication1.Optimization.LinearSystem
{
    public sealed class MINRES
    {
        #region Fields

        private const double precisionConst = 1E-25;

        public double Precision { get; private set; }

        #endregion

        #region Constructor

        public MINRES(
            double precision)
        {
            Precision = precision;
        }

        public MINRES()
            : this(precisionConst)
        { }

        #endregion

        #region Public Methods

        public Vector Solve(
            Vector[] A,
            Vector b,
            Vector startX,
            int nIter)
        {
            Vector[] symmA = A;
            Vector normb = b;
            Vector[] At = Vector.Transpose(A);

            //Symmetrize matrix
            if (!Vector.Equals(A, At))
            {
                symmA = Vector.Mult(At, A);
                normb = Vector.Mult(At, b);
            }

            Vector v0 = new Vector(b.Count());
            Vector v1 = normb - Vector.Mult(symmA, startX);
           
            double beta1 = v1.Length();
            double betaN = 0.0;
            double n = beta1;
            double c0 = 1.0;
            double c1 = 1.0;
            double s0 = 0.0;
            double s1 = 0.0;
            Vector w0 = new Vector(v1.Count());
            Vector w_1 = new Vector(v1.Count());
            Vector x = new Vector(startX);
            
            for (int i = 0; i < nIter; i++)
            {
                //Calculate Lanczos Vectors
                Vector v = (1.0 / beta1) * v1;
                Vector Av = Vector.Mult(symmA, v);
                double alpha = v * Av;
                v1 = Av - alpha * v - beta1 * v0;
                betaN = v1.Length();

                //Calculate QR factors
                double lambda = c1 * alpha - c0 * s1 * beta1;
                double p1 = Math.Sqrt(lambda * lambda + betaN * betaN);
                double p2 = s1 * alpha + c0 * c1 * beta1;
                double p3 = s0 * beta1;

                //Calculate New Givens Rotations
                c0 = c1;
                c1 = lambda / p1;

                s0 = s1;
                s1 = betaN / p1;

                //Update Solution
                Vector w = (1.0 / p1) * (v - p3 * w_1 - p2 * w0);

                x = x + c1 * n * w;
                n = -s1 * n;

                if (Math.Abs(n) < precisionConst)
                    break;

                beta1 = betaN;
                v0 = v;
                w_1 = w0;
                w0 = w;
            }

            return x;
        }

        #endregion
    }
}
