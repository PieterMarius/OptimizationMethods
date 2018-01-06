using System;

namespace ConsoleApplication1.Optimization.LinearSystem
{
    public sealed class MINRES
    {
        #region Fields

        private const double precisionConst = 1E-25;

        private double residual = double.NaN;

        public double Precision { get; private set; }
        public bool CheckSymmetry { get; private set; }

        #endregion

        #region Constructor

        public MINRES(
            double precision,
            bool checkSymmetry)
        {
            Precision = precision;
            CheckSymmetry = checkSymmetry;
        }

        public MINRES()
            : this(precisionConst, true)
        { }

        #endregion

        #region Public Methods

        public MinVector Solve(
            MinVector[] A,
            
            MinVector b,
            MinVector startX,
            int nIter)
        {
            MinVector[] symmA = A;
            MinVector normb = b;

            //Symmetrize matrix
            if (CheckSymmetry)
            {
                MinVector[] At = MinVector.Transpose(A);

                if (!MinVector.Equals(A, At))
                {
                    symmA = MinVector.Mult(At, A);
                    normb = MinVector.Mult(At, b);
                }
            }

            

            MinVector v0 = new MinVector(b.Count);
            MinVector v1 =  normb - MinVector.Mult(symmA, startX);

            double beta1 = v1.Length();
            double betaN = 0.0;
            double n = beta1;
            double c0 = 1.0;
            double c1 = 1.0;
            double s0 = 0.0;
            double s1 = 0.0;
            MinVector w0 = new MinVector(v1.Count);
            MinVector w_1 = new MinVector(v1.Count);
            MinVector x = new MinVector(startX);

            for (int i = 0; i < nIter; i++)
            {
                //Calculate Lanczos Vectors
                MinVector v = (1.0 / beta1) * v1;
                MinVector Av =  MinVector.Mult(symmA, v);
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
                MinVector w = (1.0 / p1) * (v - p3 * w_1 - p2 * w0);

                x = x + c1 * n * w;
                n = -s1 * n;

                residual = Math.Abs(n);

                if (residual < precisionConst)
                    break;

                beta1 = betaN;
                v0 = v;
                w_1 = w0;
                w0 = w;
            }

            return x;
        }

        public double GetResidual()
        {
            return residual;
        }

        #endregion
    }
}
