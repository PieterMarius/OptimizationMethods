using System;

namespace ConsoleApplication1.Optimization.LinearSystem
{
    public sealed class MINRES
    {
        #region Fields

        private const double precisionConst = 1E-30;

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

        public OptVector Solve(
            OptVector[] A,
            
            OptVector b,
            OptVector startX,
            int nIter)
        {
            OptVector[] symmA = A;
            OptVector normb = b;

            //Symmetrize matrix
            if (CheckSymmetry)
            {
                OptVector[] At = OptVector.Transpose(A);

                if (!OptVector.Equals(A, At))
                {
                    symmA = OptVector.Mult(At, A);
                    normb = OptVector.Mult(At, b);
                }
            }

            

            OptVector v0 = new OptVector(b.Count);
            OptVector v1 =  normb - OptVector.Mult(symmA, startX);

            double beta1 = v1.Length();
            double betaN = 0.0;
            double n = beta1;
            double c0 = 1.0;
            double c1 = 1.0;
            double s0 = 0.0;
            double s1 = 0.0;
            OptVector w0 = new OptVector(v1.Count);
            OptVector w_1 = new OptVector(v1.Count);
            OptVector x = new OptVector(startX);

            for (int i = 0; i < nIter; i++)
            {
                //Calculate Lanczos Vectors
                OptVector v = (1.0 / beta1) * v1;
                OptVector Av =  OptVector.Mult(symmA, v);
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
                OptVector w = (1.0 / p1) * (v - p3 * w_1 - p2 * w0);

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
